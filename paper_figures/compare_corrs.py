import pandas
import glob
import h5py
import os
import math
import numpy as np
from scipy.io import mmread
from scarlink.src.read_model import read_model
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def calc_fishers_p(r12, r13, r23, n):
    # adjust for numerical issues
    det_s = max(0, 1 - (r12**2 + r13**2 + r23**2) + 2*r12*r13*r23)

    t = (n-1)*(1+r23)
    t /= 2*((n-1)/(n-3))*det_s + (r12+r13)**2/4*math.pow(1-r23, 3)
    t = abs(r12-r13)*math.sqrt(t)
    pval = stats.t.sf(t, n-3)
    return pval

def get_gsm_gene(gsm, gene_names, rm, gene):
    # gsm: Gene score matrix with cell x gene score
    # gene_names: Names of genes in columns

    gene_gsm = np.ravel(gsm[gene_names['name'] == gene, :].todense())
    return gene_gsm[rm.test_ix]

def get_model_gene(rm, f, gene):
    # get SCARlink predictions for given gene
    w_mat = np.array(f['genes/' + gene][:])
    e = np.array([f['genes/' + gene].attrs['intercept']])
    params = [w_mat, e]
    train_alpha = float(f['genes/' + gene].attrs['alpha'])
    gex_train, gex_test = rm.get_gex_gene(gene)
    tile_gene_mat_train, tile_gene_mat_test = rm.gene_tile_matrix_scaled(gene, 
                                                    normalization_factor='ReadsInTSS')
    model_custom = rm.build_model(tile_gene_mat_test.shape[1], train_alpha)
    model_custom.set_weights(params)
    pred_vals = np.ravel(model_custom(tile_gene_mat_test.todense(), training=False).numpy())
    return pred_vals

def make_gsm_scarlink_table(dirname, out_prefix, check_hvg=True, coassay_file=""):
    # create table with Spearman correlations between
    # predicted gene expression from SCARlink and 
    # ArchR gene score matrix and observed gene expression.

    gsm_file = out_prefix + "gene_score_matrix.mtx"
    gsm_gene_info_file = out_prefix + "gene_score_matrix_gene_info.txt"
    
    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    
    if check_hvg:
        # Initially model was run on more genes. We are going to use only
        # the top 5000 highly variable genes
        # The genes are saved in hvg.txt in the scarlink_out directory
        # If file hvg.txt not found then correlation plot will be generated 
        # for all genes
        if os.path.isfile(dirname + 'hvg.txt'):
            hvg = pandas.read_csv(dirname + 'hvg.txt', header=None)[0].values.tolist()
        else:
            hvg = []

    dirname = dirname + "scarlink_out/"
    filename = out_prefix + '_values.csv'

    if os.path.isfile(filename):
        df = pandas.read_csv(filename, sep='\t')
        return df

    gsm = mmread(gsm_file).tocsr()
    gsm_gene_names = pandas.read_csv(gsm_gene_info_file, sep='\t')
    all_coef_files = glob.glob(dirname + 'coefficients*.hd5')
    # tbl = pandas.DataFrame(columns=['gene', 'model_test_corr'])
    gs_corrs = []
    model_corrs = []
    model_gs_corrs = []
    better_method = []
    gene_name = []
    method_corr = []
    for coef_file in all_coef_files:
        f = h5py.File(coef_file, mode = 'r')
        f_genes = list(f['genes/'].keys())
        if coassay_file == '':
            rm = read_model(dirname, out_file_name=coef_file.split('/')[-1], read_only=True)
        else:
            rm = read_model(dirname, out_file_name=coef_file.split('/')[-1], read_only=True, input_file_name=coassay_file)

        if check_hvg and hvg != []:
            f_genes = list(filter(lambda x: x in hvg, f_genes))

        for gene in f_genes:
            g = get_gsm_gene(gsm, gsm_gene_names, rm, gene)
            m = get_model_gene(rm, f, gene)
            _, obs = rm.get_gex_gene(gene)
            obs = np.ravel(obs.todense())
            m_obs, _ = stats.spearmanr(m, obs)
            g_obs, _ = stats.spearmanr(g, obs)
            m_g, _ = stats.spearmanr(m, g)
            gs_corrs.append(g_obs)
            model_corrs.append(m_obs)
            model_gs_corrs.append(m_g)
            fishers_corr = calc_fishers_p(m_obs, g_obs, m_g, rm.test_ix.shape[0])
            method_corr.append(fishers_corr)
            gene_name.append(gene)
            better_method.append(1 if m_obs > g_obs else -1)
        f.close()
    df = pandas.DataFrame(columns=['SCARlink corr', 'ArchR gene corr', 'gene'])
    df['SCARlink corr'] = model_corrs
    df['ArchR gene corr'] = gs_corrs
    df['SCARlink ArchR corr'] = model_gs_corrs
    df['gene'] = gene_name
    df['better_method_flag'] = better_method
    df['fishers_p'] = method_corr
    df.to_csv(filename, sep='\t', index=None)
    return df

def plot_compare_dorc_corr(dirname, out_prefix, plot_title, coassay_file='', check_hvg=True):

    df = make_dorc_scarlink_table(dirname, out_prefix, coassay_file=coassay_file, check_hvg=check_hvg)

    min_val = np.nanmin(df[['dorc gene corr', 'SCARlink corr']].values)
    max_val = np.nanmax(df[['dorc gene corr', 'SCARlink corr']].values)
    pv = df['fishers_p'].values
    pv[np.isnan(pv)] = 1
    _, pv, _, _ = multipletests(pv, method='fdr_bh')

    df['corr_significance_log'] = -np.log10(pv)
    better_method = np.ones(df.shape[0])
    better_method[(df['SCARlink corr'] < df['dorc gene corr']).values] = -1
    df['better_method_flag'] = better_method
    df['corr_significance_log'] = df['corr_significance_log']*df['better_method_flag']

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#E6A0C4", "white", "#7294D4"])
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "#E6A0C4"])
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "#7294D4"])
    fig, ax = plt.subplots(1, 3, width_ratios=[1, 0.05, 0.05], figsize=(5.5, 3.5))
    ax[0].scatter(df['dorc gene corr'], df['SCARlink corr'], c=df['corr_significance_log'], cmap=cmap, alpha=0.8, vmin=-4, vmax=4, edgecolors='white', linewidth=0.7)
    norm = plt.Normalize(-4, 4)
    norm2 = plt.Normalize(0, 4) 
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm1 = plt.cm.ScalarMappable(cmap=cmap1, norm=norm2) # "BrBG"
    sm2 = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2) # "BrBG"
    sm.set_array([])
    sm1.set_array([])
    sm2.set_array([])
    fig.colorbar(sm1, cax=ax[1], orientation='vertical')
    fig.colorbar(sm2, cax=ax[2], orientation='vertical')
    ax[0].plot([min_val-0.01, max_val+0.01], [min_val-0.01, max_val+0.01], ls='--', color='black')
    ax[0].set_xlim((min_val-0.01, max_val+0.01))
    ax[0].set_ylim((min_val-0.01, max_val+0.01))
    ax[0].set_xlabel("DORC correlation")
    ax[0].set_ylabel("SCARlink correlation")
    plt.suptitle(' '.join(plot_title.split('_')))
    plt.tight_layout()
    plt.savefig(out_prefix + '_compare_dorc.pdf', transparent=True)
    plt.show()
    plt.close()
    print("SCARlink better:", df[(df['SCARlink corr'] > df['dorc gene corr']) & (pv < 0.05)].shape[0], df[(df['SCARlink corr'] > df['dorc gene corr']) & (pv < 0.05)].shape[0]/df.shape[0], df.shape)
    print("DORC better:", df[(df['SCARlink corr'] < df['dorc gene corr']) & (pv < 0.05)].shape[0], df[(df['SCARlink corr'] < df['dorc gene corr']) & (pv < 0.05)].shape[0]/df.shape[0], df.shape)

    nan_vals = (~df['SCARlink corr'].isna()).values & (~df['dorc gene corr'].isna()).values
    _, pv = stats.wilcoxon(df['SCARlink corr'].values[nan_vals], df['dorc gene corr'].values[nan_vals], alternative='greater')
    print("Wilcoxon p-val:", pv)

def plot_compare_gsm_corr(dirname, output_prefix, plot_title, coassay_file='', check_hvg=True):

    df = make_gsm_scarlink_table(dirname, output_prefix, coassay_file=coassay_file, check_hvg=check_hvg)

    print("Number of genes:", df.shape)
    min_val = np.nanmin(df[['ArchR gene corr', 'SCARlink corr']].values)
    max_val = np.nanmax(df[['ArchR gene corr', 'SCARlink corr']].values)

    pv = df['fishers_p'].values
    pv[np.isnan(pv)] = 1
    _, pv, _, _ = multipletests(pv, method='fdr_bh')

    df['corr_significance_log'] = -np.log10(pv)
    better_method = np.ones(df.shape[0])
    better_method[(df['SCARlink corr'] < df['ArchR gene corr']).values] = -1
    df['better_method_flag'] = better_method
    df['corr_significance_log'] = df['corr_significance_log']*df['better_method_flag']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#E6A0C4", "white", "#7294D4"])
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "#E6A0C4"])
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "#7294D4"])
    fig, ax = plt.subplots(1, 3, width_ratios=[1, 0.05, 0.05], figsize=(5.5, 3.5))
    ax[0].scatter(df['ArchR gene corr'], df['SCARlink corr'], c=df['corr_significance_log'], cmap=cmap, alpha=0.8, vmin=-4, vmax=4, edgecolors='white', linewidth=0.7)
    norm = plt.Normalize(-4, 4) 
    norm2 = plt.Normalize(0, 4) 
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) # "BrBG"
    sm1 = plt.cm.ScalarMappable(cmap=cmap1, norm=norm2) # "BrBG"
    sm2 = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2) # "BrBG"
    sm.set_array([])
    sm1.set_array([])
    sm2.set_array([])
    fig.colorbar(sm1, cax=ax[1], orientation='vertical')
    fig.colorbar(sm2, cax=ax[2], orientation='vertical')
    ax[0].plot([min_val-0.01, max_val+0.01], [min_val-0.01, max_val+0.01], ls='--', color='black')
    ax[0].set_xlim((min_val-0.01, max_val+0.01))
    ax[0].set_ylim((min_val-0.01, max_val+0.01))
    ax[0].set_xlabel("ArchR GSM correlation")
    ax[0].set_ylabel("SCARlink correlation")
    plt.suptitle(' '.join(plot_title.split('_')))
    plt.tight_layout()
    plt.savefig(output_prefix + '_compare_gsm.pdf', transparent=True)
    plt.show()

    print("SCARlink better:", df[(df['SCARlink corr'] > df['ArchR gene corr']) & (pv < 0.05)].shape[0], df[(df['SCARlink corr'] > df['ArchR gene corr']) & (pv < 0.05)].shape[0]/df.shape[0], df.shape)

    print("ArchR better:", df[(df['SCARlink corr'] < df['ArchR gene corr']) & (pv < 0.05)].shape[0], df[(df['SCARlink corr'] < df['ArchR gene corr']) & (pv < 0.05)].shape[0]/df.shape[0], df.shape)
    plt.close()

    nan_vals = (~df['SCARlink corr'].isna()).values & (~df['ArchR gene corr'].isna()).values
    _, pv = stats.wilcoxon(df['SCARlink corr'].values[nan_vals], df['ArchR gene corr'].values[nan_vals], alternative='greater')
    print("Wilcoxon p-val:", pv)

def keep_dorc_genes(peak_gene):
    n_peak_gene = peak_gene.groupby('Gene').size().reset_index()
    n_peak_gene = n_peak_gene[n_peak_gene[0] >= 10]['Gene'].values
    return list(n_peak_gene)

def calc_dorc_score(peak_gene, peak_coords, peak_counts, peak_info, gene, rm, idx, dorc_genes):
    peak_gene_association = peak_gene[peak_gene['Gene'] == gene]
    if peak_gene_association.empty: return -np.ones(rm.test_ix.shape[0])

    peak_coords_gene = peak_coords[peak_coords['peak_name'].isin(peak_gene_association['peak'])]
    neg_corr = peak_gene_association[(peak_gene_association['peak'].isin(peak_coords['peak_name']))&(peak_gene_association['Corr'] < 0)]
    peak_gene_association = peak_gene[peak_gene['Gene'] == gene]
    peak_gene_association = peak_gene_association[peak_gene_association['Corr'] >= 0]
    peak_coords_gene = peak_coords[peak_coords['peak_name'].isin(peak_gene_association['peak'])]

    peak_counts_gene = np.asarray(peak_counts[peak_coords_gene.index, :].todense())[:, rm.test_ix] / np.array(peak_info.iloc[rm.test_ix]['ReadsInPeaks'])
    _, gene_test = rm.get_gex_gene(gene)
    gene_rna_norm = np.ravel(gene_test.todense())[idx]
    gene_dorc_norm = np.ravel(np.sum(peak_counts_gene, axis = 0))[idx]
    return gene_dorc_norm

def make_dorc_scarlink_table(dirname, out_prefix, check_hvg=True, coassay_file=""):
    # create table with Spearman correlations between
    # predicted gene expression from SCARlink and 
    # DORC score and observed gene expression.

    filename = out_prefix + 'dorc_values.csv'
    if os.path.isfile(filename):
        df = pandas.read_csv(filename, sep='\t')
        return df

    dirname = dirname + '/' if dirname[-1] != '/' else dirname
    if check_hvg:
        # Initially model was run on more genes. We are going to use only
        # the top 5000 highly variable genes
        # The genes are saved in hvg.txt in the scarlink_out directory
        # If file hvg.txt not found then correlation plot will be generated 
        # for all genes
        if os.path.isfile(dirname + 'hvg.txt'):
            hvg = pandas.read_csv(dirname + 'hvg.txt', header=None)[0].values.tolist()
        else:
            hvg = []

    dirname = dirname + "scarlink_out/"
    peak_gene_file = out_prefix + 'peak_gene_associations.csv'
    peak_gene = pandas.read_csv(peak_gene_file, sep=',', index_col=0)
    dorc_genes = keep_dorc_genes(peak_gene)

    peak_count_file = out_prefix + 'peak_counts.mtx'
    peak_counts = mmread(peak_count_file).tocsr()

    peak_coords_file = out_prefix + 'peak_coords.csv'
    peak_coords = pandas.read_csv(peak_coords_file, sep='\t')
    peak_coords['peak_name'] = [x['seqnames'] + '_' + str(x['start'] - 1) for i, x in peak_coords.iterrows()]

    peak_info_file = out_prefix + 'peak_info.csv'
    peak_info = pandas.read_csv(peak_info_file, sep='\t')

    all_coef_files = glob.glob(dirname + 'coefficients*.hd5')
    gs_corrs = []
    model_corrs = []
    gene_name = []
    corr_test = []
    corr_test_one_sided = []
    better_method = []
    method_corr = []
    for coef_file in all_coef_files:
        f = h5py.File(coef_file, mode = 'r')
        f_genes = list(f['genes/'].keys())

        if coassay_file == '':
            rm = read_model(dirname, out_file_name=coef_file.split('/')[-1], read_only=True)
        else:
            rm = read_model(dirname, out_file_name=coef_file.split('/')[-1], read_only=True, input_file_name=coassay_file)

        if check_hvg and hvg != []:
            f_genes = list(filter(lambda x: x in hvg, f_genes))

        for gene in f_genes:
            dorc_score = calc_dorc_score(peak_gene, peak_coords, peak_counts, peak_info, gene, rm, np.arange(rm.test_ix.shape[0]), dorc_genes)
            model_s = get_model_gene(rm, f, gene)
            _, gex_test = rm.get_gex_gene(gene)
            obs = np.ravel(gex_test.todense())
            m_obs, _ = stats.spearmanr(model_s, obs)
            if np.all(dorc_score==-1): continue
            g_obs, _ = stats.spearmanr(dorc_score, obs)
            m_g, _ = stats.spearmanr(model_s, dorc_score)

            gs_corrs.append(g_obs)
            model_corrs.append(m_obs)
            fishers_corr = calc_fishers_p(m_obs, g_obs, m_g, rm.test_ix.shape[0])
            method_corr.append(fishers_corr)
            if m_obs > g_obs: a = 'greater'
            else: a = 'lesser'
            if a == 'greater': better_method.append(1)
            else: better_method.append(-1)
            gene_name.append(gene)
        f.close()

    df = pandas.DataFrame(columns=['SCARlink corr', 'dorc gene corr', 'gene'])
    df['SCARlink corr'] = model_corrs
    df['dorc gene corr'] = gs_corrs
    df['gene'] = gene_name
    df['fishers_p'] = method_corr
    df['better_method_flag'] = better_method
    df.to_csv(filename, sep='\t', index=None)
    return df
