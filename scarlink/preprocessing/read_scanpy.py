import scanpy as sc
import os
import numpy as np
import pandas
import h5py
import tempfile
from scipy import sparse
from multiprocessing import Pool
from functools import partial

def read_chrom_sizes(filename, genome):
    chrom_sizes = dict([(x.strip().split('\t')[0], int(x.strip().split('\t')[1])) for x in  open(filename, 'r').readlines()])
    if genome == 'hg38' or genome == 'hg19':
        chrom_names = ['chr'+str(i+1) for i in range(22)] + ['chrX']
        chrom_sizes_keys = list(chrom_sizes.keys())
        for c in chrom_sizes_keys:
            if c not in chrom_names:
                chrom_sizes.pop(c)
    elif genome == 'mm10':
        chrom_names = ['chr'+str(i+1) for i in range(19)] + ['chrX']
        chrom_sizes_keys = list(chrom_sizes.keys())
        for c in chrom_sizes_keys:
            if c not in chrom_names:
                chrom_sizes.pop(c)
    elif genome == 'rn7':
        chrom_names = ['chr'+str(i+1) for i in range(20)] + ['chrX']
        chrom_sizes_keys = list(chrom_sizes.keys())
        for c in chrom_sizes_keys:
            if c not in chrom_names:
                chrom_sizes.pop(c)
    return chrom_sizes

def write_hdf5_rna(h5_file, scrna_object, selected_genes):
    scrna_object_subset = scrna_object[:, scrna_object.var_names.isin(selected_genes.gene.values)]
    df = scrna_object_subset.var.reset_index(names='gene_name')
    df.to_hdf(h5_file, key='/gene_expression/gene_info', mode='a')

    # Saving matrix information
    X = scrna_object_subset.X.copy()
    X.sort_indices()

    f = h5py.File(h5_file, mode='a')
    f.create_dataset('gene_expression/data',data= X.data)
    f.create_dataset('gene_expression/shape', data=np.flip(X.shape))
    f.create_dataset('gene_expression/indices', data=X.indices)
    f.create_dataset('gene_expression/indptr', data=X.indptr)
    f.close()
    
def write_hdf5(h5_file, scatac_object, selected_genes, gene):
    ch = selected_genes[selected_genes.gene == gene].iloc[0]['seqnames']
    start = selected_genes[selected_genes.gene == gene].iloc[0]['start']
    end = selected_genes[selected_genes.gene == gene].iloc[0]['end']
    scatac_gene = scatac_object[:, (scatac_object.var['seqnames']==ch) &
                                   (scatac_object.var['end'] >= start) &
                                   (scatac_object.var['start'] <= end)]

    # Saving matrix information.
    df = scatac_gene.var.reset_index(names='idx')[['seqnames', 'start', 'end']].copy() # .to_records(index=False)
    df.to_hdf(h5_file, key='/'+gene+'/tile_info', mode='a')
    X = scatac_gene.X.copy()
    X.sort_indices()
    f = h5py.File(h5_file, mode='a')
    f.create_dataset(gene+'/data',data= X.data)
    f.create_dataset(gene+'/shape', data=np.flip(X.shape))
    f.create_dataset(gene+'/indices', data=X.indices)
    f.create_dataset(gene+'/indptr', data=X.indptr)
    f.close()
    
# def split_write_hdf5(dirname, selected_genes, scatac_object):
def split_write_hdf5(dirname, hvg, scatac_object, max_ix, x):
    selected_genes = hvg.iloc[range(x, hvg.shape[0], max_ix)]
    with tempfile.NamedTemporaryFile(suffix=".h5", dir=dirname, mode="w+b", delete=False) as tmp_file:
        temp_h5_path = tmp_file.name

    for gene in selected_genes.gene.values:
        write_hdf5(temp_h5_path, scatac_object, selected_genes, gene)
    return temp_h5_path

def gather_tmp_h5(output_file, tmp_files):
    with h5py.File(output_file, 'w') as dest_f:
        for i, tmp_file_path in enumerate(tmp_files):
            with h5py.File(tmp_file_path, 'r') as src_f:
                object_names = list(src_f.keys())
                for obj_name in object_names:
                    src_f.copy(obj_name, dest_f, name=obj_name)

def get_gene_coordinates(gtf_file, chrom_sizes, window):
    transcripts = {}
    a = pandas.read_csv(gtf_file, sep = "\t", header = None)
    a = a[a[2].isin(['transcript', 'CDS'])]
    a = a[a[0].isin(list(chrom_sizes.keys()))]
    for i, r in a.iterrows():
        gene = r[8].split(';')[0][9:-1]
        if gene[:3] == 'MIR': continue
        if gene in transcripts:
            if r[3] < transcripts[gene][1]:
                transcripts[gene][1] = r[3]
            if r[4] > transcripts[gene][2]:
                transcripts[gene][2] = r[4]
        else:
            transcripts[gene] = [r[0], r[3], r[4]]
    transcripts = pandas.DataFrame.from_dict(transcripts).T
    transcripts.columns = ['seqnames', 'start', 'end']
    transcripts['start'] = [max(1, r['start']-window) for i,r in transcripts.iterrows()]
    transcripts['end'] = [min(chrom_sizes[r['seqnames']], r['end']+window) for i,r in transcripts.iterrows()]
    return transcripts
    
def py_write_files(scatac_out, scrna_out, out_dir, window, ncores, scale, gtf_file, chrom_sizes_file, genome):
    scatac_object = sc.read_h5ad(scatac_out)
    scrna_object = sc.read_h5ad(scrna_out)

    # Extract same cells
    scatac_object = scatac_object[scatac_object.obs_names.isin(scrna_object.obs_names)]
    scrna_object = scrna_object[scrna_object.obs_names.isin(scatac_object.obs_names)]

    # Make sure same cells are in both scrna_object and scatac_object
    scatac_object = scatac_object[np.argsort(scatac_object.obs_names)].copy()
    scrna_object = scrna_object[np.argsort(scrna_object.obs_names)].copy()

    # Normalize and scale values by  median or 10,000
    scrna_object.X = scrna_object.layers['counts'].copy()
    if scale == 'median':
        sc.pp.normalize_total(scrna_object)
    else:
        sc.pp.normalize_total(scrna_object, target_sum=1e4)

    # Remove blacklisted tiles
    if 'blacklist_tiles' in scatac_object.var.columns:
        scatac_object.X = scatac_object.X@sparse.diags((~scatac_object.var.blacklist_tiles.values).astype(int))
    if 'nFrags' not in scatac_object.obs.columns:
        scatac_object.obs['nFrags'] = scatac_object.X.sum(1)
    scatac_object.var['seqnames'] = [x.split(':')[0] for x in scatac_object.var_names]
    scatac_object.var['start'] = [int((x.split(':')[1]).split('-')[0]) for x in scatac_object.var_names]
    scatac_object.var['end'] = [int(x.split('-')[1]) for x in scatac_object.var_names]
    scatac_object.var['seqnames'] = scatac_object.var['seqnames'].astype(str)
    scatac_object.var['start'] = scatac_object.var['start'].astype(int)
    scatac_object.var['end'] = scatac_object.var['end'].astype(int)
    
    os.makedirs(out_dir, exist_ok=True)
    
    n_cells = scrna_object.shape[0]

    if n_cells == 0:
        print("No cells in common between scRNA and scATAC objects")
        exit(1)
    else:
        print(n_cells, "cells in common between scATAC and scRNA objects")
        
    chrom_sizes = read_chrom_sizes(chrom_sizes_file, genome)
    transcripts = get_gene_coordinates(gtf_file, chrom_sizes, window)
    hvg = transcripts[transcripts.index.isin(scrna_object.var_names[scrna_object.var['highly_variable']])].reset_index(names='gene')

    # create chunks
    max_ix = min(hvg.shape[0], max(hvg.shape[0]//50,50))

    tmp_files = []
    for x in range(max_ix):
        t = split_write_hdf5(out_dir, hvg, scatac_object, max_ix, x)
        tmp_files.append(t)
    # TO DO: Parallelize 
    # partial = split_write_hdf5(out_dir=out_dir, hvg=hvg, scatac_object=scatac_object, max_ix=max_ix)
    # with Pool(processes=ncores) as pool:
    #     results = pool.map(partial, range(max_ix))
    
    h5file = out_dir + '/coassay_matrix.h5'
    gather_tmp_h5(h5file, tmp_files)
    for tmp_file in tmp_files:
        os.remove(tmp_file)
        
    write_hdf5_rna(h5file, scrna_object, hvg)

    common_cols = scrna_object.obs.columns[scrna_object.obs.columns.isin(scatac_object.obs.columns)]
    if len(common_cols)==0:
        df = scrna_object.obs.merge(scatac_object.obs, left_index=True, right_index=True)
    else:
        df = scrna_object.obs.merge(scatac_object.obs.loc[:, ~scatac_object.obs.columns.isin(common_cols)], left_index=True, right_index=True)
    df.to_hdf(h5file, key='/cell_info', mode='a', format='t')

    
    # TO DO: Add LSI and KNN
    
    # Save variable genes list
    np.savetxt(out_dir + "hvg.txt", hvg['gene'].values, fmt='%s')
