import numpy as np
from scipy import stats,sparse
import math
import pandas
from sklearn.preprocessing import MinMaxScaler
from sklearn import linear_model
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
import shap

def set_gene_tile_significance_bootstrapped(x, y, w_mat, e, cell_info, celltype_col, clusters):
    """Compute standardized Shapley values for each group of cell clusters.
    
    Parameters
    ----------
    x : [[float]]
        Tile matrix for a gene.
    y : [float]
        Gene expression vector for same gene.
    w_mat : [float]
        Learned regression coefficients.
    e : float
        Learned regression bias.
    cell_info : data frame
        Cell metadata found in key cell_info in coassay_matrix.h5
    celltype_col : str
        Column in cell_info containing cell clusters.
    clusters : [str]
        List of clusters

    Returns
    -------
    zscore_d
        Dictionary of standardized z-scores for each tile and cell type
        for a given gene.
    """

    n = 500
    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    ypred = clf.predict(x)
    explainer_linear = shap.LinearExplainer(clf, x) 
    # clusters = cell_info[celltype_col].unique()
    pvals_d = {}
    shap_lst = []
    prod_lst = []
    b_clust = []
    w_mat = np.ravel(w_mat)
    shap_lst_c = []
    pv_lst = []
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        if np.sum(clust_idx) < 100: continue
        b_clust.append(clust)
        x_clust = x[clust_idx].copy()
        y_clust = y[clust_idx].copy()
        shap_lst_c = []
        ypred_subset = np.zeros((n, w_mat.shape[0]))
        ypred_minus_f = np.zeros((n, w_mat.shape[0]))
        for i in range(n):
            idx = np.random.choice(x_clust.shape[0], 50, replace=True)
            x_test_subset = x_clust[idx]
            shap_values_linear = explainer_linear.shap_values(x_test_subset)
            mean_shap = np.mean(shap_values_linear, axis=0)
            shap_lst_c.append(mean_shap)
        shap_lst.append(shap_lst_c)
    
    shap_lst = np.array(shap_lst)
    shap_lst = np.mean(shap_lst, axis=1)
    z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
    zscore_d = dict(zip(b_clust, z_score))
    for c in clusters:
        if c not in zscore_d:
            zscore_d[c] = np.zeros(w_mat.shape[0])
    return zscore_d


def set_gene_tile_significance_signed_rank(x, y, w_mat, e, cell_info, celltype_col, clusters, z_d):
    """Compute the significance of the difference in gene expression prediction when a tile
    is zero-ed out.
    
    Parameters
    ----------
    x : [[float]]
        Tile matrix for a gene.
    y : [float]
        Gene expression vector for same gene.
    w_mat : [float]
        Learned regression coefficients.
    e : float
        Learned regression bias.
    cell_info : data frame
        Cell metadata found in key cell_info in coassay_matrix.h5
    celltype_col : str
        Column in cell_info containing cell clusters.
    zscore_d : dictionary
        Dictionary of standardized z-scores for each tile and cell type
        for a given gene.
    clusters : [str]
        List of clusters

    Returns
    -------
    p_vals_d
        Dictionary of p-values for each tile and cell type
        for a given gene.
    """

    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    ypred = clf.predict(x) # .astype(np.float32)
    # clusters = cell_info[celltype_col].unique()
    b_clust = []
    p_vals_d = {}
    
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        if np.sum(clust_idx) < 100: continue
        b_clust.append(clust)
        x_clust = x[clust_idx].copy()
        ypred_clust = ypred[clust_idx].copy()
        all_idx = np.arange(x_clust.shape[1]).astype(int)
        all_p_vals = np.ones(all_idx.shape[0])
        # separately set each tile to 0 and predict                                                                                         
        for idx in all_idx:
            # check if tile accessibility is 0
            # if z_d[clust][idx] <= 0 or x_clust[:, idx].sum() == 0: continue
            if x_clust[:, idx].sum() == 0: continue
            x_clust_mod = x_clust.copy()
            x_clust_mod[:, idx] = 0
            ypred_clust_mod = clf.predict(x_clust_mod) # .astype(np.float32)
            if np.sum(ypred_clust != ypred_clust_mod) < 20: continue
            _, pv = stats.wilcoxon(ypred_clust, ypred_clust_mod, alternative='greater')
            all_p_vals[idx] = pv
        p_vals_d[clust] = all_p_vals

    for c in clusters:
        if c not in p_vals_d:
            p_vals_d[c] = np.ones(w_mat.shape[0])
    return p_vals_d
