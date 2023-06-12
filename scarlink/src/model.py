import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib import cm
import seaborn
import tensorflow as tf
import numpy as np
import pandas
import h5py
import os
import warnings
from tables import NaturalNameWarning
from scipy import stats
from scipy.sparse import csr_matrix
from sklearn.preprocessing import MaxAbsScaler
from sklearn.model_selection import RepeatedKFold, train_test_split
import tensorflow.keras.backend as K
from scarlink.src.plotExtra import plotRegion, get_fragment_counts, plot_hist, create_colormap
from scarlink.src.read_h5_and_group_cells import construct_cell_info, construct_gex_mat, get_train_test_split, get_gene_tile_matrix_group_cells
from scarlink.src.tile_significance import set_gene_tile_significance, set_gene_tile_significance_bootstrapped, set_gene_tile_significance_new_corr

warnings.filterwarnings('ignore', category=NaturalNameWarning)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

class RegressionModel:
    def __init__(self, input_file, output_dir, gtf_file = "", scatac_fragment_file = "", out_file_name = 'coefficients.hd5', group_cells = False, mode='a'):
        self.output_dir = output_dir if output_dir[-1] == '/' else output_dir + '/'
        # create output directory
        os.makedirs(self.output_dir, exist_ok = True)
        self.out_file = out_file_name 
        self.input_file_handle = pandas.HDFStore(input_file, 'r')
        # get cell info
        self.cell_info = construct_cell_info(self.input_file_handle, group_cells)        # get gene info
        self.gene_info = self.input_file_handle.select('gene_expression/gene_info')
        # construct sparse scRNA expression matrix
        self.gex_matrix = construct_gex_mat(self.input_file_handle, self.cell_info, group_cells)
        # get available list of genes
        self.gene_names = self.gene_info['gene_name']
        # get idx for training set and completely held out test set
        self.train_ix, self.test_ix = get_train_test_split(self.input_file_handle, self.gex_matrix.shape[0], random_state = 9, group_cells = group_cells)
        # create scaler object
        self.scaler = MaxAbsScaler(copy = False)
        self.alphas = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 5, 10, 20]
        # add path to GTF file
        self.gtf_file = gtf_file
        self.scatac_fragment_file = scatac_fragment_file
        self.group_cells = group_cells
        # create output file
        f = h5py.File(self.output_dir + self.out_file, mode = mode)
        if 'genes' not in f:
            g = f.create_group("genes")
            g.attrs['output_dir'] = os.path.abspath(self.output_dir)
            g.attrs['out_file'] = os.path.abspath(self.out_file)
            g.attrs['input_file_name'] = os.path.abspath(input_file)
            g.attrs['alphas'] = self.alphas
            g.attrs['gtf_file'] = os.path.abspath(self.gtf_file)
            g.attrs['scatac_fragment_file'] = os.path.abspath(self.scatac_fragment_file)
            
        f.close()
        
    def gene_tile_matrix(self, gene):
        tile_gene_mat = get_gene_tile_matrix_group_cells(self.input_file_handle, gene, self.group_cells)
        return tile_gene_mat
        
    def gene_tile_matrix_scaled_all(self, gene, normalization_factor):
        # construct sparse tile matrix for the given gene
        tile_gene_mat = self.gene_tile_matrix(gene)
        row_indices, col_indices = tile_gene_mat.nonzero()
        norm_factor = np.array(self.cell_info[normalization_factor])[row_indices]
        tile_gene_mat.data /= norm_factor 
        self.scaler.fit(tile_gene_mat)
        self.scaler.transform(tile_gene_mat)
        return tile_gene_mat
   
    def gene_tile_matrix_scaled(self, gene, normalization_factor):
        # construct sparse tile matrix for the given gene
        tile_gene_mat = self.gene_tile_matrix(gene)
        row_indices, col_indices = tile_gene_mat.nonzero()
        norm_factor = np.array(self.cell_info[normalization_factor])[row_indices]
        tile_gene_mat.data /= norm_factor

        # split into train and test
        tile_gene_mat_train = tile_gene_mat[self.train_ix]
        tile_gene_mat_test = tile_gene_mat[self.test_ix]

        # scale values
        self.scaler.fit(tile_gene_mat_train)
        self.scaler.transform(tile_gene_mat_train)
        self.scaler.transform(tile_gene_mat_test)

        return tile_gene_mat_train, tile_gene_mat_test


    def get_gex_gene(self, gene):
        # get gene expression data
        gene_gex = self.gex_matrix[:, self.gene_info['gene_name'] == gene]
        scrna_train = gene_gex[self.train_ix]
        scrna_test = gene_gex[self.test_ix]
        return scrna_train, scrna_test

    def check_gex_sparsity(self, gex_mat):
        # Check how sparse the data is 
        Y_zeros = (gex_mat.shape[0] - gex_mat.data.shape[0])/gex_mat.shape[0]
        return Y_zeros


    def find_correlation_spearman(self, model, x, y):
        y_pred = model.predict(x)
        corr, pval = stats.spearmanr(y, y_pred)
        return corr, pval

    def get_gene_tile_significance(self, gene, celltype_col):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        z_d = self.compute_gene_tile_significance(gene, celltype_col)
        dfs = []
        tiles = self.input_file_handle.select(gene + '/tile_info')
        for c in z_d:
            c_df = pandas.DataFrame(columns=['chr', 'start', 'end', celltype_col, 'z-score'])
            c_df['chr'] = tiles['seqnames']
            c_df['start'] = tiles['start'].astype(int)
            c_df['end'] = tiles['end'].astype(int)
            c_df[celltype_col] = c
            c_df['z-score'] = z_d[c]
            dfs.append(c_df)
        df = pandas.concat(dfs, axis=0)
        df['gene'] = gene
        f.close()
        return df

    def compute_gene_tile_significance(self, gene, celltype_col):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        w_mat = np.array(f['genes/' + gene][:])
        e = np.array([f['genes/' + gene].attrs['intercept']])
        train_alpha = float(f['genes/' + gene].attrs['alpha'])
        gex_train, gex_test = self.get_gex_gene(gene)
        tile_gene_mat_train, tile_gene_mat_test = self.gene_tile_matrix_scaled(gene,
                        normalization_factor='ReadsInTSS')
        x = np.array(tile_gene_mat_train.todense())
        z_d = set_gene_tile_significance_bootstrapped(x, np.ravel(gex_train.todense()), w_mat, e, self.cell_info.iloc[self.train_ix], celltype_col)
        f.close()
        return z_d 

    def build_model(self, atac_shape, a):
        # regression framework
        inputs = tf.keras.layers.Input(shape=(atac_shape,), name = 'inputA')
        out = tf.keras.layers.Dense(1, activation = tf.exp, name = 'rate', kernel_regularizer=tf.keras.regularizers.l2(a), kernel_constraint = tf.keras.constraints.NonNeg())(inputs)
        m = tf.keras.models.Model(inputs=inputs, outputs=out)
        return m

    def get_model_weight_average(self, weights):
        new_weights = list()
        for weights_list_tuple in zip(*weights): 
            new_weights.append(np.array([np.array(w).mean(axis=0) for w in zip(*weights_list_tuple)]))
        return new_weights

    
    def run_model_cross_validation(self, rna, atac, epochs, verbose, plot_loss):

        kfold = 5 
        cv = RepeatedKFold(n_splits=kfold, n_repeats = 1, random_state=9)

        best_corr = -1
        s_best_corr = -1
        best_testloss = 1000000
        flag = 0

        print("Performing " + str(int(kfold)) + "-fold cross validation...")
        if verbose: bar = progressbar.ProgressBar(max_value=len(self.alphas) * kfold)
        counter = 0
        p_aix = -1
        for aix in range(len(self.alphas)):
            a = self.alphas[aix]
            corrs = []
            s_corrs = []
            ws = []
            testlosses = []
            model_weights = []
            for train_ix, test_ix in cv.split(range(rna.shape[0])):
                rna_train = rna[train_ix]
                rna_test = rna[test_ix]
                atac_train = atac[train_ix, :]
                atac_test = atac[test_ix, :]
                opt = tf.keras.optimizers.Adam(learning_rate=0.001, epsilon = 1)
                model_custom = self.build_model(atac_train.shape[1], a)
                model_custom.compile(optimizer=opt, loss='poisson')
                hist = model_custom.fit(x=atac_train, y=rna_train, validation_data=(atac_test, rna_test), epochs=epochs, verbose=0)
                if plot_loss: plot_hist(hist, "alpha: " + str(a)) 
                s_corr, s_pval = self.find_correlation_spearman(model_custom, atac_test, rna_test)
                s_corrs.append(s_corr)
                testlosses.append(hist.history['val_loss'][-1])
                ws.append(model_custom.get_layer('rate').get_weights()[0].T.flatten())
                model_weights.append(model_custom.get_weights())
                K.clear_session()
                counter += 1
                if verbose: bar.update(counter)
                
            s_corr = np.mean(s_corrs)
            testloss = np.mean(testlosses)
            w = np.mean(ws, axis = 0)
            weights = self.get_model_weight_average(model_weights) # models)
            
            if np.isnan(s_corr):
                continue

            if s_best_corr < s_corr:
            # if testloss < best_testloss:
                best_testloss = testloss
                params = a
                s_best_corr = s_corr
                s_params = a
                best_w = weights
                flag = 1
            elif aix != 0 and p_aix == aix - 1:
                p_aix = aix
            else:
                p_aix = aix
        if flag == 0: return None, None, None
        return s_best_corr, s_params, best_w

    def test_model(self, rna, atac, w, a):
        model_custom = self.build_model(atac.shape[1], a)
        model_custom.set_weights(w)
        s_corr, s_pval = self.find_correlation_spearman(model_custom, atac, rna)
        return s_corr

    def check_if_calculated(self, gene):
        # check if regression model is already calculated
        f = h5py.File(self.output_dir + self.out_file, mode = 'a')
        cond = gene in f['genes']
        f.close()
        return cond
    
    def train_test_model(self, gene, normalization_factor = 'ReadsInTSS', max_zero_fraction = 0.9, epochs = 20, verbose = True, plot_loss = False, force = False):
        if self.check_if_calculated(gene) and not force:
            print("Gene regression model already calculated. Set force=True if you want to recalculate. Returning...")
            return
        gex_train, gex_test = self.get_gex_gene(gene)
        if self.check_gex_sparsity(gex_train) > max_zero_fraction:
            return

        tile_gene_mat_train, tile_gene_mat_test = self.gene_tile_matrix_scaled(gene, normalization_factor)
        train_corr, train_alpha, train_w = self.run_model_cross_validation(np.ravel(gex_train.todense()), tile_gene_mat_train, epochs, verbose, plot_loss)
        print("Avg. cross-validation spearman correlation", train_corr)
        print("Chosen regularization parameter:", train_alpha) 
        if train_corr is None:
            print("ERROR: Regression could not be estimated. Returning...")
            return
        test_corr = self.test_model(np.ravel(gex_test.todense()), tile_gene_mat_test, train_w, train_alpha)
        print("Spearman corr on test set:", test_corr)

        # write values
        f = h5py.File(self.output_dir + self.out_file, mode = 'a')
        if "genes/" + gene in f.keys(): del f["genes/" + gene]
        dset = f.create_dataset("genes/" + gene, data = train_w[0])
        dset.attrs['intercept'] = train_w[1][0]
        dset.attrs['spearman_correlation_train'] = train_corr
        dset.attrs['spearman_correlation_test'] = test_corr
        dset.attrs['alpha'] = train_alpha
        tiles = self.input_file_handle.select(gene + '/tile_info')
        dset.attrs['chr'] = tiles.iloc[0]['seqnames']
        dset.attrs['start'] = tiles['start'].min()
        dset.attrs['end'] = tiles['end'].max()
        dset.attrs['max_zero_fraction'] = max_zero_fraction
        dset.attrs['epochs'] = epochs
        dset.attrs['scATAC norm factor'] = normalization_factor
        f.close()


    def get_gene_corr(self, gene):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            print("ERROR: ' + gene + ' not present in " + self.output_dir + self.out_file)
            f.close()
            return
        d = {}
        d['spearman_correlation_train'] = f['genes/' + gene].attrs['spearman_correlation_train']
        d['spearman_correlation_test'] = f['genes/' + gene].attrs['spearman_correlation_test']
        return d
    
    def get_gene_coefficient(self, gene):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            print("ERROR: ' + gene + ' not present in " + self.output_dir + self.out_file)
            f.close()
            return
        w_mat = f['genes/' + gene][:]
        intercept = f['genes/' + gene].attrs['intercept']
        f.close()
        return w_mat, intercept

    def get_gene_alpha(self, gene):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            print("ERROR: ' + gene + ' not present in " + self.output_dir + self.out_file)
            f.close()
            return
        alpha = f['genes/' + gene].attrs['alpha']
        f.close()
        return alpha

    def get_gene_window_coords(self, gene):
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            print("ERROR: ' + gene + ' not present in " + self.output_dir + self.out_file)
            f.close()
            return
        chrm = f['genes/' + gene].attrs['chr']
        start = int(f['genes/' + gene].attrs['start'])
        end = int(f['genes/' + gene].attrs['end'])
        return chrm, start, end
    
    def plot_gene(self, gene, groups = 'Clusters', plot_frags = False, output_header = 'figures', to_save = False, cmap = None, save_format='png', figsize=(17, 14), sort_gex=False, show_yticks=False, plot_shap=False, shap_cmap='Blues', pvals_cmap='Blues', cluster_order=[], tilesize=500, plot_pval=False, bg_transparent=False):
        if cmap is None:
            if len(cluster_order) == 0:
                clusters = sorted(list(set(self.cell_info.dropna()[groups])))
            else:
                clusters = cluster_order
            cmap = create_colormap(clusters)
        else: 
            clusters = sorted(list(cmap.keys()))
            if len(cluster_order) != 0: clusters = cluster_order

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            print("ERROR: ' + gene + ' not present in " + self.output_dir + self.out_file)
            f.close()
            return
        w_mat = f['genes/' + gene][:] 
        chrm = f['genes/' + gene].attrs['chr']
        start = int(f['genes/' + gene].attrs['start'])
        end = int(f['genes/' + gene].attrs['end'])
        sp_corr = f['genes/' + gene].attrs['spearman_correlation_test']
        print("Spearman correlation on held out test set:", f['genes/' + gene].attrs['spearman_correlation_test'])

        if plot_pval:
            zscore_d = self.compute_gene_tile_significance(gene, groups)
            min_pvals = 0
            max_pvals = np.round(np.percentile(list(zscore_d.values()), 99)) # 6
        f.close()

        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(len(clusters) + 2, 2, width_ratios = [1, 0.3], height_ratios = [1] + [1 for i in range(len(clusters) + 1)])
        plotrange = np.arange(start + tilesize//2, end + tilesize//2, tilesize).astype(int)

        axn = plt.subplot(gs[0, 0])
        ax = []
        for i in range(len(clusters)):
            ax_a = plt.subplot(gs[i + 1, 0])
            ax_r = plt.subplot(gs[i + 1, 1])
            ax.append([ax_a, ax_r])
        ax_c = plt.subplot(gs[-1, 0])

        if plot_pval: 
            ax_shap = plt.subplot(gs[0, 1])
            s_plot = ax_shap.scatter([0, 0], [0, 0], cmap=pvals_cmap, vmin=min_pvals, vmax=max_pvals)
            ax_shap.set_xlim((1, 2))
            ax_shap.set_ylim((1, 2))
            ax_shap.axis('off')
            cb = plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min_pvals, vmax=max_pvals), cmap=pvals_cmap), ax=ax_shap, orientation='horizontal', label='z-score', location='top')
            cb.ax.xaxis.set_label_position('bottom')
            cb.ax.xaxis.set_ticks_position('bottom')
            cb.set_ticks([min_pvals, max_pvals])
            cb.set_ticklabels(['0', ">=" + str(max_pvals)]) 
        X = self.gene_tile_matrix(gene)
        Y = np.ravel(self.gex_matrix[:, self.gene_info['gene_name'] == gene].todense())
        rna_max = 0

        if sort_gex:
            means = [np.mean(Y[self.cell_info[groups] == clusters[c]]) for c in range(len(clusters))]
            c_idx = np.argsort(means)
        else:
            c_idx = list(range(len(clusters)))
            
        for c in range(len(clusters)):
            clust = clusters[c_idx[c]]
            if not plot_frags:
                X_c = np.reshape(np.sum(X[self.cell_info[groups] == clust].todense(), axis = 0), (-1, 1))
            else:
                X_c = get_fragment_counts(self.scatac_fragment_file, self.cell_info, chrm, start, end)
            X_c = X_c/X[self.cell_info[groups] == clust].shape[0]

            if cmap is None:
                ax[c][0].plot(plotrange, X_c, color = 'C' + str(c))
                seaborn.violinplot(data=Y[self.cell_info[groups] == clust], ax = ax[c][1], color = 'C' + str(c), orient='h')
            else:
                ax[c][0].plot(plotrange, X_c, color = cmap[clust])
                seaborn.violinplot(data=Y[self.cell_info[groups] == clust], ax = ax[c][1], color = cmap[clust], orient='h')
            
            if plot_pval:
                max_pvals = 4
                pvals_d = zscore_d
                ylm = ax[c][0].get_ylim()
                yrange = ylm[1] - ylm[0]
                # ordering so that higher absolute shap values are plotted last
                o_idx = np.argsort(zscore_d[clust])
                ax[c][0].scatter(x=plotrange[o_idx], y=[ylm[0] - 0.1*yrange]*zscore_d[clust].shape[0], marker='o', alpha=0.5, clip_on=False, c=pvals_d[clust][o_idx], cmap=pvals_cmap, vmin=0, vmax=max_pvals)
                ax[c][0].set_ylim(ylm)


            r_max = np.max(Y[self.cell_info[groups] == clust])
            if r_max > rna_max:
                rna_max = r_max
            ax[c][0].set_xlim((start, end)) 
            ax[c][0].set_xticks([])
            ax[c][0].set_ylabel(str(clust), rotation = 'horizontal', ha = 'right', fontsize=12)
            ax[c][1].set_yticks([])
            if c != 0: ax[c][1].spines['top'].set_visible(False) 
            ax[c][1].spines['right'].set_visible(False) 
            ax[c][1].spines['left'].set_visible(False) 
            ax[c][0].spines['top'].set_visible(False) 
            ax[c][0].spines['right'].set_visible(False) 
            ax[c][0].spines['left'].set_visible(False) 

            if not show_yticks: ax[c][0].set_yticks([])
            
            if c != len(clusters)-1: 
                ax[c][1].set_xticks([])
                ax[c][1].spines['bottom'].set_visible(False) 
                ax[c][0].spines['bottom'].set_visible(False) 


        ax_c.plot(plotrange, w_mat, color = 'black')
        ax_c.spines['top'].set_visible(False) 
        ax_c.spines['right'].set_visible(False) 

        ax_c.set_xlim((start, end))
        ax_c.set_ylabel("Coefficients", fontsize=12)

        if self.gtf_file != '':
            plotRegion(chrm, start, end, axn, self.gtf_file)
        else:
            axn.axis('off')

        plt.suptitle(gene + "(corr : " + str(round(sp_corr, 4)) + ")", fontsize = 'xx-large')
        if to_save:
            os.makedirs(self.output_dir + '/' + output_header + '/', exist_ok = True)
            plt.savefig(self.output_dir + '/' + output_header + '/' + gene + '.' + save_format, transparent=bg_transparent)
            plt.close()
            print("Saved as " + self.output_dir + '/' + output_header + '/' + gene + '.' + save_format) 
        
