import argparse
import glob
import h5py
import pandas

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output directory. Must be the same as the output directory used for scarlink_processing and scarlink.")
    parser.add_argument('-f', '--filename', required=True, type=str, help="File containing the new cell cluster labels.")
    parser.add_argument('--colname', required=True, type=str, help="Column name in input file containing the cell clusters.")
    parser.add_argument('--cellcol', required=False, type=str, help="Column name in input file with cell names. The cell names must match the cell names in Seurat and ArchR objects provided as input to scarlink_processing. If not provided, the --colname cluster is added in the order that it is provided in.")
    parser.add_argument('--force', required=False, type=bool, help="Whether to replace existing cell groupings of the same name. Default: False.")
    
    args = parser.parse_args()
    dirname = args.outdir
    if args.force is None: args.force = False
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    coassayfile = dirname + 'coassay_matrix.h5'
    df = pandas.read_csv(args.filename, sep='\t')

    # check if column present in given file
    if args.colname not in df.columns:
        print("ERROR: Invalid column name: " + args.colname)
        exit(0)
        
    f = pandas.HDFStore(coassayfile, 'r+')
    cell_info = f.select('cell_info')

    # check if column already present 
    if args.colname in cell_info.columns:
        if args.force: 
            # drop column 
            cell_info.pop(args.colname)
        else:
            print(args.colname + " already present. Set --force True to overwrite.")
            f.close()
            exit(0)

    if args.cellcol is None:
        cell_info[args.colname] = df[args.colname].values
    else:
        if args.cellcol not in df.columns:
            print("ERROR: Invalid column name: " + args.cellcol)
            f.close()
            exit(0)
        cell_info = cell_info.merge(df[[args.colname, args.cellcol]], how='left', on=args.cellcol)

    # rewrite cell_info 
    cell_info.to_csv(dirname + "cell_info.txt", sep='\t', index=None)
    f.put('cell_info', cell_info)
    f.close()

if __name__ == '__main__':
    main()
