import argparse
from scarlink.src.model import RegressionModel

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, type=str)
    parser.add_argument('-g', '--genome', required=True, type=str, choices=['mm10', 'hg38', 'hg19'])
    parser.add_argument('-p', '--proc', required=False, type=int)
    parser.add_argument('-np', '--nproc', required=False, type=int)
    parser.add_argument('-c', '--celltype', required=False, type=str)
    args = parser.parse_args()

    dirname = args.outdir
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    output_dir = dirname + 'scarlink_out/'
    input_file = dirname + 'coassay_matrix.h5'
    gtf_file = args.genome
    path = '/'.join(__file__.split('/')[:-1]) + '/data/'
    gtf_file = path + '/data/' + gtf_file + '.refGene.gtf.gz'
    if args.proc is not None:
        p_ix = args.proc-1
        if args.nproc is not None: total_procs = args.nproc
        else: total_procs = 100
    else: p_ix = None

    celltype_col = args.celltype

    rm = RegressionModel(input_file, output_dir, gtf_file=gtf_file, out_file_name = 'coefficients_' + str(p_ix) + '.hd5', group_cells = False)
    if p_ix is None:
        gene_names = rm.gene_names
    else:
        gene_names = [rm.gene_names[i] for i in range(p_ix, len(rm.gene_names), total_procs)]

    for gene in gene_names:
        rm.train_test_model(gene, normalization_factor = 'ReadsInTSS',
                            force = False, epochs = 20, verbose = False) 
        if celltype_col is not None and rm.check_if_calculated(gene):
            z = rm.compute_gene_tile_significance(gene, celltype_col)

if __name__ == '__main__':
    main()
