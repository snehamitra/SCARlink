import os
import logging
import argparse
import pandas
from datetime import datetime
from scarlink.src.model import RegressionModel

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory. Must be the same as the output directory used for scarlink_processing.")
    parser.add_argument('-g', '--genome', required=True, type=str, help="Default genome choices include mm10, hg38, hg19. Alternatively, a GTF file can be provided as input.") # , choices=['mm10', 'hg38', 'hg19']
    parser.add_argument('-p', '--proc', required=False, type=int, help="Process id. This should be set to the job index when running in parallel on a cluster.")
    parser.add_argument('-np', '--nproc', required=False, type=int, help="Number of parallel jobs. Default value is 100 when -p is provided.")
    parser.add_argument('-c', '--celltype', required=False, type=str, help="Cell type column name. This column should be present in the Seurat/ArchR object provided in scarlink_processing.")
    parser.add_argument('--sparsity', required=False, type=float, help="Maximum allowed sparsity in gene expression vector to run the regression model. Default is 0.9 meaning there can be at most 90 percent zeros in the gene expression vector.")
    parser.add_argument('--gene_list', required=False, type=str, help="Text file with gene names to run SCARlink on. Note that these genes must be a subset of highly variable genes in the Seurat object. By default SCARlink is run on all highly variable genes above the given sparsity threshold") 
    args = parser.parse_args()
    dirname = args.outdir
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    output_dir = dirname + 'scarlink_out/'
    log_dir = dirname + 'log/'
    input_file = dirname + 'coassay_matrix.h5'
    gtf_file = args.genome
    path = '/'.join(__file__.split('/')[:-1]) + '/data/'
    if gtf_file in ['mm10', 'hg19', 'hg38']:
        gtf_file = path + '/data/' + gtf_file + '.refGene.gtf.gz'
    if args.proc is not None:
        p_ix = args.proc-1
        if args.nproc is not None: total_procs = args.nproc
        else: total_procs = 100
    else: p_ix = None
    if args.sparsity is None:
        args.sparsity = 0.9
    
    celltype_col = args.celltype

    os.makedirs(log_dir, exist_ok = True)
    # create log file
    log_file_name = datetime.now().strftime('scarlink_log_' + str(p_ix) + '_%H_%M_%S_%d_%m_%Y.log')
    print("Log file: " + log_file_name)
    log_file_name = log_dir + log_file_name
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', filename=log_file_name, level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    rm = RegressionModel(input_file, output_dir, log_file_name, gtf_file=gtf_file, out_file_name = 'coefficients_' + str(p_ix) + '.hd5')

    # subset genes if gene-list is provided
    if args.gene_list is not None:
        gene_list = pandas.read_csv(args.gene_list, header=None)[0].values
        rm.gene_names = rm.gene_names[rm.gene_names.isin(gene_list)]
        print("Running SCARlink on gene set from", args.gene_list)

    if p_ix is None:
        gene_names = rm.gene_names
    else:
        gene_names = [rm.gene_names[i] for i in range(p_ix, len(rm.gene_names), total_procs)]

    for gene in gene_names:
        rm.train_test_model(gene, normalization_factor='ReadsInTSS',
                            epochs=20, verbose=False, max_zero_fraction=args.sparsity) 

    # computing Shapley values
    if celltype_col is not None:
        logging.info("Computing Shapley values for " + celltype_col + " clusters")
        for gene in gene_names:
            if rm.check_if_calculated(gene):
                z = rm.compute_gene_tile_significance(gene, celltype_col)

    rm.input_file_handle.close()

if __name__ == '__main__':
    main()
