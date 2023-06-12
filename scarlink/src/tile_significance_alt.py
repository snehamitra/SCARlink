import numpy as np
from scipy import stats,sparse
import shap
import math
import pandas
from sklearn.preprocessing import MinMaxScaler
from sklearn import linear_model
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests

def calc_pval(bootstrap_ci):
    pvals = np.ones(bootstrap_ci.confidence_interval.low.shape)
    # keep only the ones for which CI doesn't contain 0s
    # idx = np.arange(pvals.shape[0])[~((bootstrap_ci.confidence_interval.low < 0) * (bootstrap_ci.confidence_interval.high > 0))]
    idx = np.arange(pvals.shape[0])[bootstrap_ci.confidence_interval.low >= 0]

    for i in idx:
        p = stats.percentileofscore(bootstrap_ci.bootstrap_distribution[i], 0)/bootstrap_ci.bootstrap_distribution[i].shape[0]
        pvals[i] = 2*(min(p, 1-p))
    pvals[pvals==0] = 2/bootstrap_ci.bootstrap_distribution[i].shape[0]
    print(pvals)
    return pvals

def check_sparsity(x, max_sparsity=0.95):
    sparsity = 1 - np.sum(x>0, axis=0) / x.shape[0]
    non_sparse_idx = np.arange(x.shape[1])[sparsity < max_sparsity]
    return non_sparse_idx

def set_gene_tile_significance_bootstrapped(x, y, w_mat, e, cell_info, celltype_col):
    n = 500
    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    ypred = clf.predict(x)
    explainer_linear = shap.LinearExplainer(clf, x) 
    clusters = cell_info[celltype_col].unique()
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
        non_sparse_idx = check_sparsity(x_clust)
        y_clust = y[clust_idx].copy()
        shap_values_linear = explainer_linear.shap_values(x_clust)
        bootstrap_ci = stats.bootstrap((shap_values_linear[:, non_sparse_idx],), 
                                        lambda x: np.mean(x, axis=0), 
                                n_resamples=1000, confidence_level=0.95,
                         random_state=9, method='percentile')
        print(bootstrap_ci.confidence_interval.low.shape, x_clust.shape)
        pvals = np.ones(x_clust.shape[1])
        if np.sum(bootstrap_ci.confidence_interval.low >= 0) > 0:
            pvals[non_sparse_idx] = calc_pval(bootstrap_ci)

        if np.sum(pvals <= 0.1) != 0:
            # more bootstrap
            bootstrap_ci = stats.bootstrap((shap_values_linear[:, pvals <= 0.1],), 
                                        lambda x: np.mean(x, axis=0), 
                                n_resamples=50000, confidence_level=0.95,
                         random_state=9, method='percentile')
            print(bootstrap_ci.confidence_interval.low.shape, x_clust.shape)
            pvals[pvals <= 0.1] = calc_pval(bootstrap_ci)
        pvals_d[clust] = -np.log10(pvals)
    for c in clusters:
        if c not in pvals_d:
            pvals_d[c] = np.zeros(w_mat.shape[0])
    return pvals_d


