import os
import glob
import numpy as np
import h5py
import sys
import pandas
import sklearn.neighbors
from scarlink.src.model import RegressionModel 
from scarlink.src.read_model import read_model

def get_gene_gex_tiles(rm, gene):
    """Get gene expression vector and associated tile matrix.
    
    Parameters
    ----------
    rm : RegressionModel
        Regression model object.
    gene : str
        Gene for which the information is to be extracted.
    
    Returns
    -------
    gene_gex, tile_gene_mat
        Gene expression vector and tile matrix.
    """

    gene_gex = rm.gex_matrix[:, rm.gene_info['gene_name'] == gene]
    tile_gene_mat = rm.gene_tile_matrix(gene)
    row_indices, col_indices = tile_gene_mat.nonzero()
    norm_factor = np.array(rm.cell_info['ReadsInTSS'])[row_indices]
    tile_gene_mat.data /= norm_factor
    return gene_gex, tile_gene_mat

def get_y_unscaled(dirname, genes, yp_file, yo_file, smooth_vals=True, nbrs=None, all_genes=False):
    """Get predicted and observed gene expression for given set of genes.
    
    Parameters
    ----------
    dirname : str
        Output directory name. This is the output directory used as --outdir parameter
        when running scarlink.
    genes : [str]
        Set of genes for which to extract gene expression.
    yp_file : str
        Filename to save predicted gene expression.
    yo_file : str
        Filename to store observed gene expression.
    smooth_vals : bool
        Whether to smooth predicted and observed gene expression.
    nbrs : [int]
        Not implemented.
    all_genes : bool
        Whether to get predicted and observed gene expression for all genes.
    
    Returns
    -------
    y_preds_save.T, y_obs_save.T
        Predicted and observed gene expression data frames.
    """

    all_coef_files = glob.glob(dirname + 'coefficients*.hd5')
    if os.path.isfile(yp_file) and os.path.isfile(yo_file):
        y_pred_save = pandas.read_csv(yp_file, sep='\t')
        y_obs_save = pandas.read_csv(yo_file, sep='\t') 
        return y_pred_save, y_obs_save
    y_preds_save = []
    y_obs_save = []
    gene_order = []
    all_ws = []
    for coef_file in all_coef_files:
        f = h5py.File(coef_file, mode = 'r')
        f_genes = list(f['genes/'].keys())
        f.close()
        rm = read_model(dirname, out_file_name = coef_file.split('/')[-1])
        flag = 0 
        for gene in f_genes:
            if not all_genes and gene not in genes: continue
            corrs = rm.get_gene_corr(gene)
            gene_gex, tile_gene_mat = get_gene_gex_tiles(rm, gene)

            w, e = rm.get_gene_coefficient(gene)
            y_pred = np.ravel(np.exp(np.dot(tile_gene_mat.todense(), w) + e))
            corrs = rm.get_gene_corr(gene)

            gene_gex = np.ravel(gene_gex.todense())

            gene_gex = np.ravel(gene_gex)
            y_pred = np.ravel(y_pred)

            if smooth_vals:
                y_pred = np.mean(np.take(y_pred, nbrs), axis=1)
                gene_gex = np.mean(np.take(gene_gex, nbrs), axis=1)

            y_preds_save.append(y_pred)
            y_obs_save.append(gene_gex)
            gene_order.append(gene)

    y_preds_save = pandas.DataFrame(y_preds_save, index=gene_order)
    y_obs_save = pandas.DataFrame(y_obs_save, index=gene_order)
    y_preds_save.T.to_csv(yp_file, sep='\t', index=None)
    y_obs_save.T.to_csv(yo_file, sep='\t', index=None)

    return y_preds_save.T, y_obs_save.T

def smooth_vals(out_dir, lsi, k):
    """Smooth predicted and observed gene expression over k nearest neighbor graph.
    
    Parameters
    ----------
    out_dir : str
        Directory in which SCARlink regression outputs are going to be saved. The function creates 
        directory <output dir>/scarlikn_out in which results are saved.
    lsi : matrix
        LSI matrix.
    k : int
        k-nearest neighbors for kNN graph.
    
    Returns
    -------
    yp_new, yo_new
        Predicted and observed gene expression data frames.
    """

    u_pred_file = out_dir + '/pred_unsmooth.csv'
    u_obs_file = out_dir + '/obs_unsmooth.csv'
    yp_df, yo_df = get_y_unscaled(out_dir, [], u_pred_file,
                                  u_obs_file, all_genes=True, 
                                  smooth_vals=False)
    yp = yp_df.values
    yo = yo_df.values
    yp_new = np.zeros(yp.shape).astype(np.float32)
    yo_new = np.zeros(yo.shape).astype(np.float32)
    
    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=k, n_jobs=-1, metric='cosine')
    nn.fit(lsi)
    dists, neighs = nn.kneighbors(lsi)
    for i in range(lsi.shape[0]):
        yp_new[i] = np.mean(yp[neighs[i]], axis=0)
        yo_new[i] = np.mean(yo[neighs[i]], axis=0)
    
    yp_new = pandas.DataFrame(yp_new, columns=yp_df.columns)
    yo_new = pandas.DataFrame(yo_new, columns=yo_df.columns)
    return yp_new, yo_new
