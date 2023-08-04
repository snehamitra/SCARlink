import argparse
from rpy2.robjects.packages import STAP

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scrna', required=True, type=str)
    parser.add_argument('--scatac', required=True, type=str)
    parser.add_argument('-o', '--outdir', required=True, type=str)
    parser.add_argument('--window', required=False, type=int)
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
    rfunc.write_files(archr_out, seurat_out, out_dir, args.window)
    
if __name__ == '__main__':
    main()
