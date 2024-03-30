import argparse
from rpy2.robjects.packages import STAP
from configparser import ConfigParser

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scrna', required=True, type=str, help="Seurat object of scRNA-seq data. The cell names must match the cell names in the scATAC-seq object.")
    parser.add_argument('--scatac', required=True, type=str, help="ArchR object of scATAC-seq data. The cell names must match the cell names in scRNA-seq object.")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory. The same output directory must be used later to run scarlink and scarlink_tiles.")
    parser.add_argument('--window', required=False, type=int, help="Number of bases to consider beyond the gene body both upstream and downstream. Default is 250000.")
    parser.add_argument('-nc', '--ncores', required=False, type=int, help="Number of cores to parallelize the preprocessing step. Default is 5.")
    parser.add_argument('--scale', required=False, type=str, choices=['median', '10k', 'none'], help="scRNA-seq normalization scaling factor. Default is none; meaning re-normalization will not be done. The normalized matrix in the Seurat object will be used.")
    args = parser.parse_args()

    rscript_dir = '/'.join(__file__.split('/')[:-1]) + '/'
    with open(rscript_dir+'read_seurat_archr.R', 'r') as f:
        string = f.read()
    rfunc = STAP(string, "rfunc")
    
    archr_out = args.scatac
    seurat_out = args.scrna
    out_dir = args.outdir
    out_dir = out_dir + '/' if out_dir[-1] != '/' else out_dir
    if args.window is None: 
        args.window = 250000
    if args.ncores is None:
        args.ncores = 5
    if args.scale is None:
        args.scale = 'none'
    print("Using " + str(args.ncores) + " cores")

    # create config 
    config_object = ConfigParser()
    config_object["PREPROCESSINFO"] = {
        "scrna": seurat_out,
        "scatac": archr_out, 
        "outdir": out_dir, 
        "window": args.window, 
        "ncores": args.ncores,
        "scale": args.scale
    }


    rfunc.write_files(archr_out, seurat_out, out_dir, args.window, args.ncores, args.scale)

    # write config file
    with open(out_dir + "config.ini", 'w') as conf:
        config_object.write(conf)
    
if __name__ == '__main__':
    main()
