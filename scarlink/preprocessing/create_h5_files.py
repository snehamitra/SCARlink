import argparse
from configparser import ConfigParser
from scarlink.preprocessing.read_scanpy import py_write_files

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scrna', required=True, type=str, help="Seurat/AnnData object of scRNA-seq data. The cell names must match the cell names in the scATAC-seq object. If AnnData then the scATAC-seq object must also be AnnData.")
    parser.add_argument('--scatac', required=True, type=str, help="ArchR/AnnData object of scATAC-seq data. The cell names must match the cell names in scRNA-seq object. If AnnData then the scRNA-seq object must also be AnnData.")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory. The same output directory must be used later to run scarlink and scarlink_tiles.")
    parser.add_argument('-g', '--genome', required=True, type=str, help="Default genome choices include mm10, mm39, hg38, hg19. Alternatively, a GTF file can be provided as input.") # , choices=['mm10', 'mm39', 'hg38', 'hg19']
    parser.add_argument('--chrom-sizes', required=False, type=str, help="chrom.sizes file. Required when using custom genome.")
    parser.add_argument('--norm-factor-scatac', required=False, type=str, help="Normalization factor for tile matrix. Default is nFrags. Replace with other factor. Original version used ReadsInTSS from ArchR.")
    parser.add_argument('--window', required=False, type=int, help="Number of bases to consider beyond the gene body both upstream and downstream. Default is 250000.")
    parser.add_argument('-nc', '--ncores', required=False, type=int, help="Number of cores to parallelize the preprocessing step. Default is 5.")
    parser.add_argument('--scale', required=False, type=str, choices=['median', '10k'], help="scRNA-seq normalization scaling factor. Default is 10k")
    args = parser.parse_args()

    scatac_out = args.scatac
    scrna_out = args.scrna
    out_dir = args.outdir
    out_dir = out_dir + '/' if out_dir[-1] != '/' else out_dir
    gtf_file = args.genome
    genome = args.genome
    norm_factor = 'nFrags' if args.norm_factor_scatac is None else args.norm_factor_scatac
    path = '/'.join(__file__.split('/')[:-1]) + '/../data/'
    if gtf_file in ['mm10', 'hg19', 'hg38', 'mm39']:
        chrom_sizes = path + gtf_file + '.chrom.sizes'
        gtf_file = path + gtf_file + '.refGene.gtf.gz'
    else:
        chrom_sizes = args.chrom_sizes
        genome = 'custom'
    if args.window is None: 
        args.window = 250000
    if args.ncores is None:
        args.ncores = 5
    if args.scale is None:
        args.scale = '10k'
    print("Using " + str(args.ncores) + " cores")
    if scrna_out[-4:] == 'h5ad' and scatac_out[-4:] == 'h5ad':
        input_format = 'anndata'
    else:
        input_format = 'rds'
    # create config 
    config_object = ConfigParser()
    config_object["PREPROCESSINFO"] = {
        "scrna": scrna_out,
        "scatac": scatac_out, 
        "outdir": out_dir, 
        "window": args.window, 
        "ncores": args.ncores,
        "scale": args.scale,
        "gtf_file": gtf_file,
        "chrom_sizes": chrom_sizes,
        "input_format": input_format,
        "norm_factor_scatac": norm_factor
    }


    if scrna_out[-4:] == 'h5ad' and scatac_out[-4:] == 'h5ad':
        print("Input files in anndata format")
        py_write_files(scatac_out, scrna_out, out_dir, args.window, args.ncores, args.scale, gtf_file, chrom_sizes, genome)
    else:
        from rpy2.robjects.packages import STAP
        print("Input files in Seurat/ArchR format")
        rscript_dir = '/'.join(__file__.split('/')[:-1]) + '/'
        with open(rscript_dir+'read_seurat_archr.R', 'r') as f:
            string = f.read()
        rfunc = STAP(string, "rfunc")
        rfunc.write_files(scatac_out, scrna_out, out_dir, args.window, args.ncores, args.scale)

    # write config file
    with open(out_dir + "config.ini", 'w') as conf:
        config_object.write(conf)
    
if __name__ == '__main__':
    main()
