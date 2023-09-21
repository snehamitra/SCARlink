import argparse
import glob
import h5py
import pandas
from statsmodels.stats.multitest import multipletests
from scarlink.src.read_model import read_model

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory. Must be the same as the output directory used for scarlink_processing and scarlink.")
    parser.add_argument('-c', '--celltype', required=True, type=str, help="Cell type column name. This column should be present in the Seurat/ArchR object provided in scarlink_processing and scarlink must have been run with the same cell type column name.")
    args = parser.parse_args()

    dirname = args.outdir
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    output_dir = dirname + 'scarlink_out/'
    celltype_col = args.celltype

    all_coef_files = sorted(glob.glob(output_dir + 'coefficients*.hd5'))
    
    dfs = []
    for coef_file in all_coef_files:
        f = h5py.File(coef_file, mode = 'r')
        f_genes = list(f['genes/'].keys())
        rm = read_model(output_dir, out_file_name=coef_file.split('/')[-1],
                        read_only=True)
        for gene in f_genes:
            df = rm.get_gene_tile_significance(gene, celltype_col)
            dfs.append(df)
    df_all = pandas.concat(dfs)
    df_all['start'] = df_all['start'].astype(int)
    df_all['end'] = df_all['end'].astype(int)
    _, adj_pval, _, _ = multipletests(df_all['p-value'].values, method='fdr_bh')
    df_all['FDR'] = adj_pval
    df_all.to_csv(dirname + 'gene_linked_tiles_' + celltype_col + '.csv.gz', sep='\t', index=None)
    rm.input_file_handle.close()
    print("Saved output:", dirname + 'gene_linked_tiles_' + celltype_col + '.csv.gz')
    
if __name__ == '__main__':
    main()
