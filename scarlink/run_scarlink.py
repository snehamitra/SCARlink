import sys
from scarlink.src.model import RegressionModel

def main():
    if len(sys.argv) < 3:
        print("Usage: scarlink <out dir> <genome> <segment --optional>")
        print("Example: scarlink outDir hg38")
        exit(0)

    dirname = (sys.argv)[1] 
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    output_dir = dirname + 'scarlink_out/'
    input_file = dirname + 'coassay_matrix.h5'
    gtf_file = (sys.argv)[2]
    path = '/'.join(__file__.split('/')[:-1]) + '/data/'
    gtf_file = path + '/data/' + gtf_file + '.refGene.gtf.gz'
    if len(sys.argv) >= 5:
        slurm_ix = int((sys.argv)[3])-1
        slurm_total_procs = int((sys.argv)[4])
    else: slurm_ix = None

    celltype_col = 'celltype'

    rm = RegressionModel(input_file, output_dir, gtf_file=gtf_file, out_file_name = 'coefficients_' + str(slurm_ix) + '.hd5', group_cells = False)
    if slurm_ix is None:
        gene_names = rm.gene_names
    else:
        gene_names = [rm.gene_names[i] for i in range(slurm_ix, len(rm.gene_names), slurm_total_procs)]

    for gene in gene_names:
        rm.train_test_model(gene, normalization_factor = 'ReadsInTSS', force = False, epochs = 20, verbose = False) 
        if celltype_col is not None and rm.check_if_calculated(gene):
            z = rm.compute_gene_tile_significance(gene, celltype_col)
            
if __name__ == '__main__':
    main()
