from rpy2.robjects.packages import STAP
import sys

def main():
    rscript_dir = '/'.join(__file__.split('/')[:-1]) + '/'
    with open(rscript_dir+'read_seurat_archr.R', 'r') as f:
        string = f.read()
    rfunc = STAP(string, "rfunc")
    
    archr_out = (sys.argv)[1]
    seurat_out = (sys.argv)[2]
    out_dir = (sys.argv)[3]

    rfunc.write_files(archr_out, seurat_out, out_dir)
    
if __name__ == '__main__':
    main()
