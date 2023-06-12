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

def bootstraping(x, y, clf, cell_info, celltype_col):
    n = 500
    clusters = sorted(list(set(cell_info.dropna()[celltype_col])))
    explainer_linear = shap.LinearExplainer(clf, x)
    pvals_d = {}
    b_clust = []
    shap_lst = []
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        if np.sum(clust_idx) < 100: continue
        b_clust.append(clust)
        x_clust = x[clust_idx].copy()
        y_clust = y[clust_idx].copy()
        shap_lst_c = []
        for i in range(n):
            idx = np.random.choice(x_clust.shape[0], 50, replace=False)
            x_train_subset = np.delete(x_clust, idx, axis=0)
            x_test_subset = x_clust[idx]
            x_train_subset_mean = x_train_subset.mean(0)
            x_test_subset_mean = x_test_subset.mean(0)
            #explainer_linear = shap.LinearExplainer(clf, x_clust) # x_train_subset)
            shap_values_linear = explainer_linear.shap_values(x_test_subset)
            mean_shap = np.mean(shap_values_linear, axis=0)
            shap_lst_c.append(mean_shap)
        shap_lst.append(shap_lst_c)
    shap_lst = np.array(shap_lst)
    for clust_idx in range(len(b_clust)):
        clust = b_clust[clust_idx]
        pvals_d[clust] = np.ones(shap_lst.shape[2])
        for idx in range(pvals_d[clust].shape[0]):
            if np.all(shap_lst[clust_idx][:, idx] <= 0): continue
            row_mean = np.mean(np.delete(shap_lst[clust_idx], idx, axis=1))
            col_mean = np.mean(np.delete(shap_lst[:, :, idx], clust_idx, axis=0))
            _, p_row = stats.ttest_1samp(shap_lst[clust_idx][:, idx], popmean=max(0, row_mean), alternative='greater')
            _, p_col = stats.ttest_1samp(shap_lst[clust_idx][:, idx], popmean=max(0, col_mean), alternative='greater')
            # _, p = stats.combine_pvalues([p_row, p_col], method='pearson')
            pvals_d[clust][idx] = p_row
    return pvals_d

def set_gene_tile_significance_new_corr(x, y, w_mat, e, cell_info, celltype_col):
    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    clusters = cell_info[celltype_col].unique()
    
    explainer_linear = shap.LinearExplainer(clf, x) # _clust) # x_train_subset)
    shap_values_linear = explainer_linear.shap_values(x)
    z_scores = stats.zscore(np.ravel(shap_values_linear)).reshape(shap_values_linear.shape)
    z_all = np.mean(z_scores, axis=0)
    c_all = np.array([stats.spearmanr(np.ravel(x[:, ix]), y)[0] for ix in range(x.shape[1])])
    pvals_d = {}
    p_lst = []
    bg_lst = []
    b_clust = []
    corr_d = {}
    z_d = {}
    z_celltype_d = {}
    corr_all_d = {}
    z_all_d = {}
    z_lst = []
    # generate background
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        x_clust = x[clust_idx]
        c = np.array([stats.spearmanr(np.ravel(x_clust[:, ix]), y[clust_idx])[0] for ix in range(x_clust.shape[1])])
        corr_d[clust] = c
        z_d[clust] = np.mean(z_scores[clust_idx], axis=0)
        corr_all_d[clust] = c_all
        z_all_d[clust] = z_all
        z_celltype_d[clust] = stats.zscore(np.ravel(np.mean(shap_values_linear[clust_idx], axis=0)))
        z_lst.append(list(np.mean(shap_values_linear[clust_idx], axis=0)))
    z_lst = np.array(z_lst)
    z_lst = stats.zscore(np.ravel(z_lst)).reshape(z_lst.shape)
    z_lst_d = dict(zip(clusters, z_lst))
    return corr_d, z_d, corr_all_d, z_all_d, z_celltype_d, z_lst_d

def set_gene_tile_significance(x, y, w_mat, e, cell_info, celltype_col):
    clf = linear_model.PoissonRegressor()
    clf.fit(x, y)
    clf.coef_ = np.ravel(w_mat)
    clf.intercept_ = e
    clusters = cell_info[celltype_col].unique()
    
    explainer_linear = shap.LinearExplainer(clf, x) # _clust) # x_train_subset)
    pvals_d = {}
    print(x.shape)
    p_lst = []
    bg_lst = []
    b_clust = []
    
    # generate background
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        x_clust = x[clust_idx]
        idx = np.random.choice(x_clust.shape[0], 100, replace=True)
        shap_values_linear = explainer_linear.shap_values(x_clust)
        bg_lst.extend(list(shap_values_linear))
    bg_lst = np.ravel(np.array(bg_lst))
    
    # estimate significance
    for clust in clusters:
        print(clust)
        clust_idx = (cell_info[celltype_col] == clust).values
        x_clust = x[clust_idx]
        shap_values_linear = explainer_linear.shap_values(x_clust)
        pv_clust = []
        pvals = np.ones(w_mat.shape[0])
        shap_0 = np.mean(shap_values_linear, axis=0)
        pvals = np.ones(w_mat.shape[0])
        shap_0 = np.mean(shap_values_linear, axis=0)
        shap_idx = np.arange(w_mat.shape[0])[shap_0 > 0]
        for s_idx in shap_idx:
            _, pv = stats.ttest_ind(np.ravel(shap_values_linear[:, s_idx]), 
                                   bg_lst) # , alternative='greater')
            pvals[s_idx] = pv
        pvals_d[clust] = pvals
    print(pvals_d)
    return pvals_d

def set_gene_tile_significance_old(x, y, w_mat, e, cell_info, celltype_col):
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
    w_mat = np.ravel(w_mat)
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        x_clust = x[clust_idx].copy()
        ypred_clust = ypred[clust_idx].copy()
        shap_values_linear = explainer_linear.shap_values(x_clust)
        mean_shap = np.mean(shap_values_linear, axis=0)
        shap_lst.append(mean_shap)
        
    shap_lst = np.array(shap_lst)
    z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
    zscore_d = dict(zip(clusters, z_score))
    #pvals_d = bootstraping(x, y, clf, cell_info, celltype_col)
    print(list(pvals_d.keys()))
    print(list(zscore_d.keys()))
    return pvals_d, zscore_d

# def set_gene_tile_significance(x, y, w_mat, e, cell_info, celltype_col):
#     dist = getattr(stats, 'expon')
#     w_mat = np.ravel(w_mat)
#     params = dist.fit(w_mat)
#     arg = params[:-2]
#     loc = params[-2]
#     scale = params[-1]
    
#     cdf_w = stats.expon.cdf(np.ravel(w_mat), loc=loc, scale=scale)
#     x_df = pandas.DataFrame(x)
#     x_df[celltype_col] = cell_info[celltype_col].values
#     x_df = x_df.groupby(celltype_col).mean()

#     scaler = MinMaxScaler()
#     x_scaled = scaler.fit_transform(x_df)
#     w_x_scaled = x_scaled * cdf_w
#     w_celltype_d = dict(zip(list(x_df.index), list(w_x_scaled)))
#     w_all_d = dict(zip(list(x_df.index), list(cdf_w)*len(x_df.index)))
#     return w_all_d, w_celltype_d

def compute_mad(scores):
    scores_median = np.median(scores)
    scores_median_absolute_deviation = np.abs(scores - scores_median)
    median_median_absolute_deviations = np.median(scores_median_absolute_deviation)
    scores_mad = median_median_absolute_deviations * 1.4826
    scores_deviation = (scores - scores_median) / scores_mad
    return scores

def set_gene_tile_significance_bootstrapped_trylater(x, y, w_mat, e, cell_info, celltype_col):
    n = 10000
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
    ttest_lst = []
    for clust in clusters:
        clust_idx = (cell_info[celltype_col] == clust).values
        if np.sum(clust_idx) < 100: continue
        b_clust.append(clust)
        #print(clust_idx.shape, x.shape)
        x_clust = x[clust_idx].copy()
        y_clust = y[clust_idx].copy()
        #shap_lst_c = []
        #for i in range(n):
        ypred_subset = np.zeros((n, w_mat.shape[0]))
        ypred_minus_f = np.zeros((n, w_mat.shape[0]))
        ##########################################################
        # for i in range(n):
        #     idx = np.random.choice(x_clust.shape[0], 50, replace=True)
        #     x_test_subset = x_clust[idx]
        #     shap_values_linear = explainer_linear.shap_values(x_test_subset)
        #     mean_shap = np.mean(shap_values_linear, axis=0)
        #     shap_lst_c.append(mean_shap)
        ###########################################################
        idx = np.random.choice(x_clust.shape[0], n, replace=True)
        x_test_subset = x_clust[idx]
        shap_values_linear = explainer_linear.shap_values(x_test_subset)
        shap_lst_c.append(shap_values_linear)
    shap_lst_c = np.array(shap_lst_c)
    z_lst = stats.zscore(np.ravel(shap_lst_c)).reshape(shap_lst_c.shape)
    print("z shape:", z_lst.shape)
    for i in range(z_lst.shape[0]):
        ttest = np.ones(z_lst.shape[2])
        for j in range(z_lst.shape[2]):
            ttest[j] = 1 - np.sum(z_lst[i, :, j] >= 2)/z_lst.shape[1]
        ttest_lst.append(ttest)
    ttest_lst = np.array(ttest_lst)
    ttest_lst = multipletests(np.ravel(ttest_lst), method='fdr_bh')[1]
    print(ttest_lst)
    print((ttest_lst < 0.05).sum())
    exit(0)
    # zscore_d = dict(zip(b_clust, z_score))
    ttest_d = dict(zip(b_clust, ttest_lst))
    for c in clusters:
        if c not in ttest_d:
            # zscore_d[c] = np.zeros(w_mat.shape[0])
            ttest_d[c] = np.zeros(w_mat.shape[0])
    #return zscore_d, ttest_d
    return ttest_d, ttest_d

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
        #print(clust_idx.shape, x.shape)
        x_clust = x[clust_idx].copy()
        y_clust = y[clust_idx].copy()
        shap_lst_c = []
        #for i in range(n):
        ypred_subset = np.zeros((n, w_mat.shape[0]))
        ypred_minus_f = np.zeros((n, w_mat.shape[0]))
        ##########################################################
        for i in range(n):
            idx = np.random.choice(x_clust.shape[0], 50, replace=True)
            x_test_subset = x_clust[idx]
            shap_values_linear = explainer_linear.shap_values(x_test_subset)
            mean_shap = np.mean(shap_values_linear, axis=0)
            shap_lst_c.append(mean_shap)
        ###########################################################
        # idx = np.random.choice(x_clust.shape[0], n, replace=True)
        # x_test_subset = x_clust[idx]
        # shap_values_linear = explainer_linear.shap_values(x_test_subset)
        # ttest = np.ones(shap_values_linear.shape[1])
        # for i in range(ttest.shape[0]):
        #     # ttest[i] = min(np.sum(shap_values_linear[:, i] > 0), np.sum(shap_values_linear[:, i] < 0))/shap_values_linear.shape[0]*2
        #     ttest[i] = np.sum(shap_values_linear[:, i] > 0)/shap_values_linear.shape[0]*2
        # ###########################################################
        shap_lst.append(shap_lst_c)
    
    shap_lst = np.array(shap_lst)
    # # # ttest_lst = np.ones((shap_lst.shape[0], shap_lst.shape[2]))
    # # for i in range(ttest_lst.shape[0]):
    # #     for j in range(ttest_lst.shape[1]):
    # #         if np.mean(shap_lst[i, :, j]) <= 0: continue
    # #         # ttest_lst[i][j] = stats.ttest_1samp(np.ravel(shap_lst[i, :, j])[1], popmean=0)[1]
    # #         # ttest_lst[i][j] = ztest(np.ravel(shap_lst[i, :, j])[1], value=0)
    # #         ttest_lst[i][j] = min(np.sum(shap_lst[i, :, j] > 0), np.sum(shap_lst[i, :, j] < 0))/shap_lst.shape[1]*2
    # # print(ttest_lst)
    # # #ttest_lst = [[stats.ttest_1samp(np.ravel(shap_lst[c, :, t]), popmean=0)[1] for t in range(shap_lst.shape[2])] for c in range(shap_lst.shape[0])]
    # # ttest_lst = np.array(ttest_lst)
    

    # z_score = np.mean(stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape), axis=1)
    shap_lst = np.mean(shap_lst, axis=1)
    # # shap_mad = compute_mad(np.ravel(shap_lst)).reshape(shap_lst.shape)
    # # shap_mad_d = dict(zip(b_clust, shap_mad))
    z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
    # z_score = np.array([stats.zscore(x) for x in shap_lst])
    zscore_d = dict(zip(b_clust, z_score))
    # ttest_d = dict(zip(b_clust, ttest_lst))
    for c in clusters:
        if c not in zscore_d:
            zscore_d[c] = np.zeros(w_mat.shape[0])
            # ttest_d[c] = np.zeros(w_mat.shape[0])
    return zscore_d


# def set_gene_tile_significance_bootstrapped_old(x, y, w_mat, e, cell_info, celltype_col):
#     n = 500
#     clf = linear_model.PoissonRegressor()
#     clf.fit(x, y)
#     clf.coef_ = np.ravel(w_mat)
#     clf.intercept_ = e
#     ypred = clf.predict(x)
#     explainer_linear = shap.LinearExplainer(clf, x) 
#     clusters = cell_info[celltype_col].unique()
#     pvals_d = {}
#     shap_lst = []
#     prod_lst = []
#     b_clust = []
#     w_mat = np.ravel(w_mat)
#     shap_lst_c = []
#     pv_lst = []
#     for clust in clusters:
#         clust_idx = (cell_info[celltype_col] == clust).values
#         if np.sum(clust_idx) < 100: continue
#         b_clust.append(clust)
#         #print(clust_idx.shape, x.shape)
#         x_clust = x[clust_idx].copy()
#         y_clust = y[clust_idx].copy()
#         shap_lst_c = []
#         #for i in range(n):
#         ypred_subset = np.zeros((n, w_mat.shape[0]))
#         ypred_minus_f = np.zeros((n, w_mat.shape[0]))
#         for i in range(n):
#             idx = np.random.choice(x_clust.shape[0], 50, replace=False)
#             x_test_subset = x_clust[idx]
#             shap_values_linear = explainer_linear.shap_values(x_test_subset)
#             mean_shap = np.mean(shap_values_linear, axis=0)
#             shap_lst_c.append(mean_shap)

#         #     ## compare without feature
#         #     ypred_subset[i] = np.mean(np.array([clf.predict(x_test_subset)]*w_mat.shape[0]), axis=1)
#         #     for j in range(w_mat.shape[0]):
#         #         x_test_subset_j = x_test_subset.copy()
#         #         x_test_subset_j[:, j] = 0
#         #         ypred_minus_f[i][j] = np.mean(clf.predict(x_test_subset_j))
#         #     ypred_minus_f[i] = np.array([np.mean(np.exp(np.log(clf.predict(x_test_subset)) - w_mat[j]*x_test_subset[:, j])) for j in range(w_mat.shape[0])])
#         #     # v = np.log(clf.predict(x_test_subset))[:, None] - x_test_subset * w_mat
#         #     # ypred_minus_f[i] = np.mean(np.exp(v), axis=0)
#         # neq_ix = np.sum(ypred_subset != ypred_minus_f, axis=0)
#         # # get at least 100 nonequal samples
#         # neq_ix = np.arange(w_mat.shape[0])[neq_ix > 100]
#         # print("neq:", neq_ix.shape, neq_ix)

#         # pv_celltype = np.ones(w_mat.shape[0])
#         # _, pv = stats.wilcoxon(ypred_subset[:, neq_ix], ypred_minus_f[:, neq_ix], alternative='greater', axis=0)
#         # pv_celltype[neq_ix] = pv

#         # pv_lst.append(list(pv_celltype))
#         shap_lst.append(shap_lst_c)

#     shap_lst = np.array(shap_lst)
#     # pv_lst = np.array(pv_lst)
#     # z_score = np.mean(stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape), axis=1)
#     # zscore_d = dict(zip(b_clust, z_score))

#     # rank_shap = stats.rankdata(shap_lst).reshape(shap_lst.shape)
#     # areas = {}
#     # print(shap_lst.shape, rank_shap.shape)
#     # pb = np.ravel(rank_shap)
#     # B = pb.shape[0]
#     # for clust_idx in range(len(b_clust)):
#     #     clust = b_clust[clust_idx]
#     #     areas[clust] = np.zeros(x.shape[1])
#     #     for idx in range(x.shape[1]):
#     #         pf = np.ravel(rank_shap[clust_idx][:, idx])
#     #         F = pf.shape[0]
#     #         areas[clust][idx] = (pf.sum()/F - pb.sum()/B)/(B+F)

#     shap_lst = np.mean(shap_lst, axis=1)
#     # # shap_mad = compute_mad(np.ravel(shap_lst)).reshape(shap_lst.shape)
#     # # shap_mad_d = dict(zip(b_clust, shap_mad))
#     z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
#     # z_score = np.array([stats.zscore(x) for x in shap_lst])
#     zscore_d = dict(zip(b_clust, z_score))
#     z_score_all = stats.zscore(np.ravel(np.mean(shap_lst, axis=0)))
#     zscore_all_d = dict(zip(b_clust, [z_score_all]*len(b_clust)))
#     # pvals_d = dict(zip(b_clust, pv_lst))
#     for c in clusters:
#         if c not in zscore_d:
#             zscore_d[c] = np.zeros(w_mat.shape[0])
#             # areas[c] = np.zeros(w_mat.shape[0])
#             # shap_mad_d[c] = np.zeros(w_mat.shape[0])
#             zscore_all_d[c] = np.zeros(w_mat.shape[0])
#             # pvals_d[c] = np.ones(w_mat.shape[0])
#     return zscore_d, zscore_all_d # , pvals_d # shap_mad_d


# def set_gene_tile_significance(x, y, w_mat, e, cell_info, celltype_col):
#     clf = linear_model.PoissonRegressor()
#     clf.fit(x, y)
#     clf.coef_ = np.ravel(w_mat)
#     clf.intercept_ = e
    
#     explainer_linear = shap.LinearExplainer(clf, x) 
#     clusters = cell_info[celltype_col].unique()
#     pvals_d = {}
#     shap_lst = []
#     prod_lst = []
#     w_mat = np.ravel(w_mat)
#     for clust in clusters:
#         clust_idx = (cell_info[celltype_col] == clust).values
#         x_clust = x[clust_idx].copy()
#         y_clust = y[clust_idx].copy()
#         # idx = np.random.choice(x_clust.shape[0], 50, replace=False)
#         # x_train_subset = np.delete(x_clust, idx, axis=0)
#         # x_test_subset = x_clust
#         # x_train_subset_mean = x_train_subset.mean(0)
#         # x_test_subset_mean = x_test_subset.mean(0)
#         shap_values_linear = explainer_linear.shap_values(x_clust)
#         mean_shap = np.mean(shap_values_linear, axis=0)
#         # prod_shap = np.multiply(w_mat, np.mean(x_clust, axis=0))

#         prod_shap = np.zeros(w_mat.shape[0])
#         x_clust_mean = np.mean(x_clust, axis=0)
#         orig = np.exp(np.multiply(w_mat, x_clust_mean) + e).sum()
#         for j in range(w_mat.shape[0]):
#             x_minus_j = x_clust_mean.copy()
#             x_minus_j[j] = 0
#             y_minus_j = np.exp(np.multiply(w_mat, x_minus_j) + e).sum()
#             prod_shap[j] = orig - y_minus_j
#         shap_lst.append(mean_shap)
#         prod_lst.append(prod_shap)
#     shap_lst = np.array(shap_lst)
#     prod_lst = np.array(prod_lst)

#     # z_score = stats.zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)
#     z_score = stats.zscore(np.ravel(prod_lst)).reshape(prod_lst.shape)
#     pvals = stats.norm.sf(np.ravel(z_score)).reshape(shap_lst.shape)
#     pvals[z_score < 0] = 1
#     pvals_d = dict(zip(clusters, pvals))
#     zscore_d = dict(zip(clusters, z_score))
#     return pvals_d, zscore_d

