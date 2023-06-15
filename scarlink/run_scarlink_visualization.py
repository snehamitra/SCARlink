import os
import pandas
import argparse
import scarlink.src.visualization as scv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, type=str)
    parser.add_argument('--genes', required=True, type=str)
    parser.add_argument('-c', '--celltype', required=True, type=str)
    args = parser.parse_args()
    
    # check if genes are in file
    if os.path.isfile(args.genes):
        genes = pandas.read_csv(args.genes, header=None).values.tolist()
        genes = sum(genes, [])
    else:
        genes = args.genes.split(',')

    scarlink_out = scv.get_scarlink_output(args.outdir)
    scv.plot_scarlink_output(scarlink_out, genes=genes,
                             celltype=args.celltype)
    
if __name__ == '__main__':
    main()
