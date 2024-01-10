import os
os.environ['KMP_WARNINGS'] = 'off'
import matplotlib.pyplot as plt
import numpy as np
import pandas 
import anndata as ad
import scanpy as sc
import sklearn.preprocessing
from scipy import stats
import sklearn.neighbors
from sklearn.cluster import AgglomerativeClustering
from scarlink.src.get_smoothed_pred_obs import smooth_vals
import warnings
warnings.filterwarnings("ignore", message="Transforming to str index.")
warnings.filterwarnings("ignore", message="An input array is constant; the correlation coefficient is not defined.") 
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )

def create_object(dirname, smooth_k=10, use_hvg=True, celltype_col='celltype', umap=None, lsi_file=''):
    """Create SCARlink object for estimating chromatin potential. The 
    function creates two AnnData objects for predicted and observed 
    gene expression.
    
    Parameters
    ----------
    dirname : str
        Output directory name. This is the output directory used as --outdir parameter
        when running scarlink.
    smooth_k : int
        Number of nearest neighbors to smooth gene expression over.
    use_hvg : bool
        Used in development mode only.
    celltype_col : str
        The column in the cell_info data frame inside coassay_matrix.h5
        containing the cell groupings.
    umap : matrix
        User-provided UMAP. If not provided, a force directed layout is 
        computed.
    lsi_file : str
        Path to file containing LSI of scATAC-seq. By default SCARlink
        saves the lsi matrix during scarlink_preprocessing. However, an 
        alternate LSI matrix can be used.
    
    Returns
    -------
    d
        Dictionary with predicted and observed gene expression.
    """

    cell_info_file = dirname + '/cell_info.txt'
    if lsi_file=='': lsi_file = dirname + '/scatac_LSI.csv'
    out_dir = dirname + '/scarlink_out/'
    lsi = pandas.read_csv(lsi_file, sep='\t').values
    yp, yo = smooth_vals(out_dir, lsi, smooth_k)
    cell_info = pandas.read_csv(cell_info_file, sep="\t")

    yp.index = yp.index.astype(str)
    yo.index = yo.index.astype(str)

    adata_pred = ad.AnnData(yp)
    adata_obs = ad.AnnData(yo)

    adata_pred.obs = cell_info
    adata_obs.obs = cell_info
    
    ##
    if not use_hvg:
        print("Using all genes to estimate chromatin potential")
    else:
        hvg = pandas.read_csv(dirname + '/hvg.txt', sep='\t', header=None)[0].values
        adata_obs = adata_obs[:, adata_obs.var_names.isin(hvg)].copy()
        adata_pred = adata_pred[:, adata_pred.var_names.isin(hvg)].copy()
    ##

    d = {}
    d['dirname'] = dirname
    d['pred'] = adata_pred
    d['obs'] = adata_obs
    d['cell_info'] = cell_info
    d['celltype_col'] = celltype_col
    
    # use user provided embedding
    if umap is not None:
        d['umap'] = umap 
        d['obs'].obsm['X_umap'] = d['umap']
    adata_scv = d['obs'].copy()
    adata_scv.obsm['X_pca'] = lsi
    if umap is not None: adata_scv.obsm['X_umap'] = d['umap']
    sc.pp.neighbors(adata_scv, n_neighbors=20, n_pcs=lsi.shape[1])

    if umap is None:
        sc.tl.paga(adata_scv, groups=celltype_col)
        sc.pl.paga(adata_scv)
        sc.tl.draw_graph(adata_scv, layout='fa', init_pos='paga', iterations=5000)
        d['umap_key'] = 'draw_graph_fa'
    else:
        d['umap_key'] = 'umap'
    d['obs'] = adata_scv
    return d

# clustering function
def cluster_genes(d, n_clust=2):
    """Cluster genes expression matrix using hierarchical clustering.
    Note that here the genes are clustered and not the cells.
    
    Parameters
    ----------
    d : dictionary
        Dictionary with predicted and observed gene expression.
    n_clust : int
        Number of clusters.
    
    Returns
    -------
    d
        Dictionary with gene clusters.
    """

    yo_file = d['dirname'] + '/scarlink_out/obs_unsmooth.csv'
    yo = pandas.read_csv(yo_file, sep="\t").T
    yo = yo.iloc[yo.index.isin(d['pred'].var_names), ]

    cell_info_file = d['dirname'] + '/cell_info.txt'
    cell_info = pandas.read_csv(cell_info_file, sep="\t")
    adata = ad.AnnData(yo.T)
    adata.obs.index = adata.obs.index.astype(str)
    adata.obs = cell_info
    
    model = AgglomerativeClustering(n_clusters=n_clust, metric='cosine', linkage='complete') 
    yhat = model.fit_predict(sklearn.preprocessing.MinMaxScaler().fit_transform(adata.X).T)
    adata.var['hierarchical'] = pandas.Categorical([str(int(x)) for x in yhat])
    
    adata_o = d['obs'][:, np.argsort(adata.var['hierarchical'])]
    adata_o = adata_o[np.argsort(adata_o.obs[d['celltype_col']]), :].copy()
    sc.pp.normalize_total(adata_o)
    sc.pp.scale(adata_o)
    uniq_c = adata.var['hierarchical'].unique()
    var_names = {}
    for c in uniq_c:
        var_names[c] = adata.var_names[adata.var['hierarchical'] == c].tolist()
        
    sc.pl.heatmap(adata_o, var_names=var_names, 
                  groupby=d['celltype_col'], cmap='plasma', show_gene_labels=False,
                  figsize=(3, 10), standard_scale='obs')
    return var_names

def filter_genes(d, batch, genes):
    """Extract a the gene expression for a subset of provided 
    genes from the predict and observed AnnData objects.
    
    Parameters
    ----------
    d : dictionary
        Dictionary with predicted and observed gene expression.
    batch : str
        Not implemented.
    genes : [str]
        List of genes to extract.
    
    Returns
    -------
    d
        Dictionary with predicted and observed gene expression
        for a subset of genes.
    """

    if batch is None:
        batch_key = ''
    else:
        batch_key = '_batch_' + str(batch)
    d['pred' + batch_key].obs.index = d['pred' + batch_key].obs.index.astype(str)
    d['obs' + batch_key].obs.index = d['obs' + batch_key].obs.index.astype(str)
    d['pred' + batch_key] = d['pred' + batch_key][:, d['pred' + batch_key].var_names.isin(genes)]
    d['obs' + batch_key] = d['obs' + batch_key][:, d['obs' + batch_key].var_names.isin(genes)]
    return d

def smooth_arrows(x, y, u, v, smooth_w=50, min_count=5, take_top=0):
    """Smooth arrows in each window.
    
    Parameters
    ----------
    x : [float]
        Starting x-coordinates for the arrows.
    y : [float]
        Starting y-coordinates for the arrows.
    u : [float]
        Length of arrows in x-axis.
    v : [float]
        Length of arrows in y-axis.
    smooth_w : int
        Create smooth_w X smooth_w windows to smooth arrows over.
    min_count : int
        Minimum number of arrows in a window required to smooth.
    take_top : int
        Number of top arrows to select per window.
    
    Returns
    -------
    x_ms, y_ms, Ex_ms, Ey_ms
        Smoothed arrow coordinates and arrow sizes.
    """

    x_ms = []
    y_ms = []
    Ex_ms = []
    Ey_ms = []

    x_window = (np.max(x) - np.min(x)) / smooth_w
    y_window = (np.max(y) - np.min(y)) / smooth_w

    for i in np.arange(np.min(x), np.max(x), x_window):
        for j in np.arange(np.min(y), np.max(y), y_window):
            idx = (x > i) & (x < (i + x_window)) & (y > j) & (y < (j + y_window))
            if np.sum(idx) < min_count:
                continue
            if take_top == 0:
                x_m = np.mean(x[idx])
                y_m = np.mean(y[idx])
                Ex_m = np.mean(u[idx])
                Ey_m = np.mean(v[idx])

            else:
                top_idx = np.argsort(w[idx])[-take_top:]
                x_m = np.mean(x[idx][top_idx])
                y_m = np.mean(y[idx][top_idx])
                Ex_m = np.mean(u[idx][top_idx])
                Ey_m = np.mean(v[idx][top_idx])
            x_ms.append(x_m)
            y_ms.append(y_m)
            Ex_ms.append(Ex_m)
            Ey_ms.append(Ey_m)

    x_ms = np.array(x_ms)
    y_ms = np.array(y_ms)
    Ex_ms = np.array(Ex_ms)
    Ey_ms = np.array(Ey_ms)
    return x_ms, y_ms, Ex_ms, Ey_ms

def calc_velocity(d, pred_key='y_pred_scaled_filtered', 
                     obs_key='y_obs_scaled_filtered', max_per_cell=10, 
                     umap_key='umap', metric='correlation', batch=None):
    """Calculate chromatin potential.
    
    Parameters
    ----------
    d : dictionary
        Dictionary with predicted and observed gene expression.
    pred_key : str
        Key in dictionary containing predicted gene expression.
    obs_key : str
        Key in dictionary containing observed gene expression.
    max_per_cell : int
        Number of neighbors to compute chromatin potential over.
    umap_key : str
        Key to access UMAP or FDL coordinates.
    metric : str
        Metric to use to estimate similarity between predicted and 
        observed gene expression.
    batch : str
        Not implemented.
    
    Returns
    -------
    umap, V, M
        Arrow coordinates and length of arrows.
    """

    yp = d[pred_key].values 
    yo = d[obs_key].values
    umap = d['obs'].obsm['X_' + umap_key] \
                    if batch is None \
                    else d['scvelo_batch_' + str(batch)].obsm['X_' + umap_key]

    
    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=max_per_cell, n_jobs=-1, metric=metric) 
    nn.fit(yp)
    dists, neighs = nn.kneighbors(yo)
    us = []
    vs = []
    cutoff = 0.95
    umap_new_unsmooth = np.zeros(umap.shape)
    umap_new = np.zeros(umap.shape)
    arrow_length = np.zeros(umap.shape[0])
    arrow_length_unsmooth = np.zeros(umap.shape[0])
    arrow_length_new = np.zeros(arrow_length.shape)
    
    M = np.zeros((yp.shape[0], yo.shape[0]))
    M[:] = np.inf

    for gx in range(yp.shape[0]):
        umap_new_unsmooth[gx] = umap[neighs[gx]].mean(0)
        M[gx][neighs[gx]] = dists[gx]
        
    for gx in range(yp.shape[0]):
        arrow_length_unsmooth[gx] = ((umap[gx][0] - umap_new[gx][0])**2 + (umap[gx][1] - umap_new[gx][1])**2)**0.5
    umap_new_unsmooth[arrow_length_unsmooth > np.quantile(arrow_length_unsmooth, cutoff)] = umap[arrow_length_unsmooth > np.quantile(arrow_length_unsmooth, cutoff)]


    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=15, n_jobs=-1, metric='euclidean')
    nn.fit(umap)
    dists, neighs = nn.kneighbors(umap)

    for gx in range(umap.shape[0]):
        umap_new[gx] = umap_new_unsmooth[neighs[gx]].mean(0)
        
    for gx in range(yp.shape[0]):
        arrow_length[gx] = ((umap[gx][0] - umap_new[gx][0])**2 + (umap[gx][1] - umap_new[gx][1])**2)**0.5

    umap_new[arrow_length > np.quantile(arrow_length, cutoff)] = umap[arrow_length > np.quantile(arrow_length, cutoff)]
    dX = umap_new - umap
    
    V = dX 
    return(umap, V, M)

def get_corrs(d, pred_key='pred', obs_key='obs'):
    """Compute correlation between predicted and observed
    gene expression for each gene.
    
    Parameters
    ----------
    d : dictionary
        Dictionary with predicted and observed gene expression.
    pred_key : str
        Dictionary key for predicted gene expression.
    obs_key : str
        Dictionary key for observed gene expression.
    
    Returns
    -------
    corrs
        Correlations for each gene.
    """

    scaler = sklearn.preprocessing.StandardScaler()
    yp = scaler.fit_transform(d[pred_key].X)
    yo = scaler.fit_transform(d[obs_key].X)
    corrs = list(map(lambda x: stats.pearsonr(yp[:, x], yo[:, x])[0], np.arange(yp.shape[1])))
    corrs = np.array(corrs)
    return corrs

def chrom_pot(d_orig, batch=None, umap_key='umap', max_per_cell=10,
              metric='cosine', scale_max_abs_val=10,
             gene_corr_cutoff=0, scaling='standard'):
    """Compute chromatin potential.
    
    Parameters
    ----------
    d_orig : dictionary
        Dictionary with predicted and observed gene expression.
    batch : str
        Not implemented.
    umap_key : str
        Key to access UMAP or FDL coordinates.
    max_per_cell : int
        Number of neighbors to compoute chromatin potential over.
    metric : str
        Metric to use to estimate similarity between predicted and 
        observed gene expression.
    scale_max_abs_val : float
        Trimming scaled values.
    gene_corr_cutoff : float
        Subset genes based on minimum correlation between predicted and observed
        gene expression to consider for chromatin potential.
    scaling : str
        Can be 'standard' or 'minmax'
    
    Returns
    -------
    x, v, d, M
        Return arrows, arrow sizes, and dictionary with predicted and 
        observed gene expression.
    """

    d = d_orig.copy()
    if scaling == 'standard':
        scaler = sklearn.preprocessing.StandardScaler()
    elif scaling == 'minmax':
        scaler = sklearn.preprocessing.MinMaxScaler()
    if batch == None: 
        pred_key = 'pred'
        obs_key = 'obs'
    else:
        pred_key = 'pred_batch_' + str(batch)
        obs_key = 'obs_batch_' + str(batch)
    yp = scaler.fit_transform(d[pred_key].X)
    yo = scaler.fit_transform(d[obs_key].X)

    tmp = pandas.DataFrame(yp, columns=d[pred_key].var_names)
    yp[yp > scale_max_abs_val] = scale_max_abs_val
    yp[yp < -scale_max_abs_val] = -scale_max_abs_val
    yo[yo > scale_max_abs_val] = scale_max_abs_val
    yo[yo < -scale_max_abs_val] = -scale_max_abs_val

    corrs = get_corrs(d, pred_key, obs_key)
    corrs_rows = []
    for i in range(yp.shape[0]):
        idx = np.argsort(yo[i])
        idx = idx[corrs[idx] > gene_corr_cutoff]
        c, _ = stats.pearsonr(yp[i][idx], yo[i][idx])
        corrs_rows.append(c)

    corrs_rows = np.array(corrs_rows)
    if batch is None: cell_info_key = 'cell_info'
    else: cell_info_key = 'cell_info_batch_' + str(batch)

    celltypes = d[cell_info_key][d['celltype_col']].unique()

    idx = np.arange(yp.shape[1])
    idx = idx[corrs[idx] > gene_corr_cutoff]
    # print("Percentage of genes lost:", (yp.shape[1]-idx.shape[0])/yp.shape[1]*100)
    if batch is None:
        pred_filtered_key = 'y_pred_scaled_filtered'
        obs_filtered_key = 'y_obs_scaled_filtered'
    else:
        pred_filtered_key = 'y_pred_scaled_filtered_batch_' + str(batch)
        obs_filtered_key = 'y_obs_scaled_filtered_batch_' + str(batch)

    d[pred_filtered_key] = pandas.DataFrame(yp[:, idx], columns=d['pred'].var_names[idx])
    d[obs_filtered_key] = pandas.DataFrame(yo[:, idx], columns=d['obs'].var_names[idx])
    x, v, M = calc_velocity(d, pred_key=pred_filtered_key, 
                     obs_key=obs_filtered_key, umap_key=umap_key, batch=batch, max_per_cell=max_per_cell, 
                              metric=metric)

    return x, v, d, M

def plot_arrows(d_data, genes=[], smooth_w=40, figsize=(10,8), max_per_cell=10, gene_corr_cutoff=0, n_arrows=500, min_count=5, **kwargs):
    """Plot chromatin potential arrows.
    
    Parameters
    ----------
    d_data : dictionary
        Dictionary with predicted and observed gene expression.
    genes : [str]
        Gene subset to use for estimating chromatin potential.
    smooth_w : int
        Create smooth_w X smooth_w windows to smooth arrows over.
    figsize : (float, float)
        Size of output plot.
    max_per_cell : int
        Number of neighbors to compoute chromatin potential over.
    gene_corr_cutoff : float
        Subset genes based on minimum correlation between predicted and observed
        gene expression to consider for chromatin potential.
    n_arrows : int
        Number of arrows to plot.
    min_count : int 
        Minimum number of arrows in a window required to smooth.
    **kwargs 
        Additional arguments passed to quiver().
    """

    batch = None
    umap_key = d_data['umap_key']
    d = d_data.copy()
    if genes != []:
        d = filter_genes(d, batch, genes)
    d[umap_key] = d_data['obs'].obsm['X_' + umap_key]
    x, v, d, M = chrom_pot(d, batch=batch, umap_key=umap_key, max_per_cell=max_per_cell, 
                       metric='correlation', gene_corr_cutoff=gene_corr_cutoff,
                           scaling='minmax')
    r_idx = np.random.choice(np.arange(x.shape[0]),
                             size=n_arrows, replace=False) if x.shape[0] >= n_arrows \
                             else np.arange(x.shape[0])
    E_norm = np.sqrt(v[:, 0]**2 + v[:, 1]**2)[r_idx]

    scv_key = 'obs'
    fig = plt.figure(figsize=figsize)
    if umap_key == 'draw_graph_fa':
        ax = sc.pl.draw_graph(d[scv_key], color=d['celltype_col'], show=False, ax=plt.gca(), frameon=False)
    else:
        ax = sc.pl.umap(d[scv_key], color=d['celltype_col'], show=False, ax=plt.gca(), frameon=False)

    x_smooth, y_smooth, u_smooth, v_smooth = smooth_arrows(x[:, 0], x[:, 1], 
                                                           v[:, 0], v[:, 1], 
                                                           smooth_w=smooth_w, min_count=min_count)
    E_norm_smooth = np.sqrt(u_smooth**2 + v_smooth**2)
    plt.quiver(x_smooth, y_smooth,
               u_smooth/E_norm_smooth, v_smooth/E_norm_smooth, 
               angles='xy', scale_units='xy', **kwargs)
