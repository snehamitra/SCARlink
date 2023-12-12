import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib import cm
import logging
import seaborn
import tensorflow as tf
import numpy as np
import pandas
import h5py
import sys
import os
import warnings
from tables import NaturalNameWarning
from scipy import stats
from scipy.sparse import csr_matrix
from sklearn.preprocessing import MaxAbsScaler
from sklearn.model_selection import RepeatedKFold, train_test_split
import tensorflow.keras.backend as K
from scarlink.src.plotExtra import plotRegion, get_fragment_counts, plot_hist, create_colormap
from scarlink.src.read_h5_and_group_cells import construct_cell_info, construct_gex_mat, get_train_test_split, get_gene_tile_matrix_group_cells, write_significance, read_sparse_significance
from scarlink.src.tile_significance import set_gene_tile_significance_bootstrapped, set_gene_tile_significance_signed_rank

warnings.filterwarnings('ignore', category=NaturalNameWarning)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

class RegressionModel:
    """SCARlink regression model to predict gene expression from tile-based accessibility
    matrix using regularized Poisson regression.

    Parameters
    ----------
    input_file : str
        coassay_matrix.h5 file generated by scarlink_preprocessing. This file name is pulled 
        from the output directory provided by the --outdir parameter. The output directory 
        contains file coassay_matrix.h5 generated by scarlink_preprocessing.
    output_dir : str
        Directory in which SCARlink regression outputs are going to be saved. The function creates 
        directory <output dir>/scarlink_out in which results are saved.
    log_file_name : str
        Name of log file. Log files are stored in <output dir>/log/.
    gtf_file : str
        GTF file for the genome used. SCARlink provides genome annotations for hg38, hg19, and mm10.
        If a different genome is used, then this parameter would include the path to the GTF file
        of that genome.
    scatac_fragment_file : str
        Path to the fragment file for plotting accessibility. This is not completely implemented.
        Right now the output plots include tile based counts and not counts for individual bases
        that can be plotted using fragment file. Plotting tile counts are faster than plotting 
        fragment counts.
    out_file_name : str
        Name of the regression output file. The output file is saved in <output dir>/scarlink_out/. 
        The output file name has the format coefficients_<process_number>.hd5. The <process_id> is None
        if SCARlink is run sequentially on all genes in coassay_matrix.h5. If SCARlink is submitted as
        a parallel job on the cluster then process_number would be the index for a given split.
    group_cells : bool
        Not implemented. Flag to check whether to group cells into pseuo-bulk counts instead of 
        running the model on single-cell data.
    mode : str
        I/O mode for out_file_name. Default is append mode.

    Example
    -------
    >>> input_file = dirname + 'coassay_matrix.h5'
    >>> output_dir = dirname + 'scarlink_out/'
    >>> log_dir = dirname + 'log/'
    >>> log_file_name = datetime.now().strftime('scarlink_log_' + str(p_ix) + '_%H_%M_%S_%d_%m_%Y.log')
    >>> log_file_name = log_dir + log_file_name
    >>> logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', filename=log_file_name, level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    >>> p_ix = None
    >>> rm = RegressionModel(input_file, output_dir, log_file_name, gtf_file=gtf_file, out_file_name = 'coefficients_' + str(p_ix) + '.hd5')

    Notes
    -----
    See further details on https://github.com/snehamitra/SCARlink
    """

    def __init__(self, input_file, output_dir, log_file_name='', gtf_file = "", scatac_fragment_file = "", out_file_name = 'coefficients.hd5', group_cells = False, mode='a'):
        
        self.output_dir = output_dir if output_dir[-1] == '/' else output_dir + '/'
        # create output directory
        os.makedirs(self.output_dir, exist_ok = True)
        self.out_file = out_file_name 
        self.log_file_name = log_file_name
            
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
        """Extract tile-based accessibility matrix centered on given gene.
        The accessibility matrix is not normalized.
        
        Parameters
        ----------
        gene : str
            Gene name for which tile-matrix is to be extracted.
        group_cells : bool
            Not implemented. Flag to check whether to group cells into pseuo-bulk counts instead of 
            running the model on single-cell data.
        
        Returns
        -------
        tile_gene_mat
            Sparse tile accessibility matrix centered on gene body.
        """

        tile_gene_mat = get_gene_tile_matrix_group_cells(self.input_file_handle, gene, self.group_cells)
        return tile_gene_mat
        
    def gene_tile_matrix_scaled_all(self, gene, normalization_factor):
        """Extract tile-based accessibility matrix centered on given gene.
        The accessibility matrix is normalized using the normalization factor.

        Parameters
        ----------
        gene : str
            Gene name for which tile-matrix is to be extracted.
        normalization_factor : str
            Column name of normalization factor present in the key
            cell_info inside coassay_matrix.h5.
        
        Returns
        -------
        tile_gene_mat
            Sparse tile accessibility matrix centered on gene body.
        """

        # construct sparse tile matrix for the given gene
        tile_gene_mat = self.gene_tile_matrix(gene)
        row_indices, col_indices = tile_gene_mat.nonzero()
        norm_factor = np.array(self.cell_info[normalization_factor])[row_indices]
        tile_gene_mat.data /= norm_factor 
        self.scaler.fit(tile_gene_mat)
        self.scaler.transform(tile_gene_mat)
        return tile_gene_mat
   
    def gene_tile_matrix_scaled(self, gene, normalization_factor):
        """Extract tile-based accessibility matrix centered on given gene.
        The accessibility matrix is normalized using the normalization factor
        and min-max-scaled. The tile matrix is then split into train and test
        set based on train and test indices.
        
        Parameters
        ----------
        gene : str
            Gene name for which tile-matrix is to be extracted.
        normalization_factor : str
            Column name of normalization factor present in the key
            cell_info inside coassay_matrix.h5.
        
        Returns
        -------
        tile_gene_mat_train, tile_gene_mat_test
            Sparse tile accessibility matrix centered on gene body
            split by train and test set.
        """

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
        """Extract gene expression vecotr for given gene from coassay_matrix.h5. 
        The gene expression vector is then split into train and test set based 
        on train and test indices.
        
        Parameters
        ----------
        gene : str
            Gene name for which gene expression vector is to be extracted.
        
        Returns
        -------
        scrna_train, scrna_test
            Gene expression vector split by train and test set.
        """

        # get gene expression data
        gene_gex = self.gex_matrix[:, self.gene_info['gene_name'] == gene]
        scrna_train = gene_gex[self.train_ix]
        scrna_test = gene_gex[self.test_ix]
        return scrna_train, scrna_test

    def check_gex_sparsity(self, gex_mat):
        """Compute fraction of zeros in given gene expression vector.
        
        Parameters
        ----------
        gex_mat : sparse vector
            Gene expression vector.
        
        Returns
        -------
        Y_zeros
            Fraction of zeros.
        """
        # Check how sparse the data is 
        Y_zeros = (gex_mat.shape[0] - gex_mat.data.shape[0])/gex_mat.shape[0]
        return Y_zeros


    def find_correlation_spearman(self, model, x, y):
        """Compute Spearman correlation between predicted gene expression 
        and observed gene expression
        
        Parameters
        ----------
        model : TensorFlow model
            Gene-level model.
        x : matrix
            Tile matrix centered on gene body.
        y : vector
            Observed gene expression vector.

        Returns
        -------
        corr, pval
            Spearman correlation between model(x) and y and associated
            p-value.
        """
        y_pred = model.predict(x, verbose=None)
        corr, pval = stats.spearmanr(y, y_pred)
        return corr, pval

    def get_gene_tile_significance(self, gene, celltype_col):
        """For each tile in the tile matrix centered around the gene,
        compute Shapley values for each cell cluster in held-out 
        test set and assess the significance of the difference between 
        predicted gene expression with and without a given tile in any
        given cell cluster.
        
        Parameters
        ----------
        gene : str
            Gene to run the analysis on.
        celltype_col : str
            The column in the cell_info data frame inside coassay_matrix.h5
            containing the cell groupings.

        Returns
        -------
        df
            Data frame where each row is a cell cluster and a unique
            tile coordinate and its associated z-score and p-value.
        """

        clusters = sorted(self.cell_info[celltype_col].unique().tolist())
        k = 'tile_significance/' + celltype_col + '/' + gene
        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        sp_corr = f['genes/' + gene].attrs['spearman_correlation_test']
        w = f['genes/' + gene][:]
        m_z = read_sparse_significance(f, k, 'z-score')
        m_p = read_sparse_significance(f, k, 'p-value')

        df_z = pandas.DataFrame(m_z.todense(), columns=clusters)
        df_p = pandas.DataFrame(10**(-np.array(m_p.todense())), columns=clusters)
        tiles = self.input_file_handle.select(gene + '/tile_info')
        tiles = tiles[['seqnames', 'start', 'end']].rename(columns={'seqnames': 'chr'})
        _, tile_acc_test = self.gene_tile_matrix_scaled(gene, 'ReadsInTSS')
        tile_acc_sum = list(np.ravel(np.mean(tile_acc_test, axis=0)))
        tiles['test_acc'] = tile_acc_sum
        tiles['test_acc_sparsity'] = list(np.ravel(np.sum((tile_acc_test > 0), axis=0) / tile_acc_test.shape[0]))

        df_z = pandas.concat([df_z, tiles], axis=1)
        df_z['regression_coef'] = w
        df_z_long = pandas.melt(df_z, id_vars=['chr', 'start', 'end', 'test_acc', 'test_acc_sparsity', 'regression_coef'], value_vars=clusters,
                                var_name=celltype_col, value_name='z-score')
        df_p = pandas.concat([df_p, tiles], axis=1)
        df_p_long = pandas.melt(df_p, id_vars=['chr', 'start', 'end', 'test_acc', 'test_acc_sparsity'], value_vars=clusters,
                                var_name=celltype_col, value_name='p-value')

        df = df_z_long.merge(df_p_long, on=['chr', 'start', 'end', 'test_acc', 'test_acc_sparsity']+[celltype_col])
        df['gene'] = gene
        df['Spearman corr'] = sp_corr
        return df
    
    def compute_gene_tile_significance(self, gene, celltype_col):
        """For each tile in the tile matrix centered around the gene,
        compute Shapley values for each cell cluster in held-out 
        test set and assess the significance of the difference between 
        predicted gene expression with and without a given tile in any
        given cell cluster.
        
        Parameters
        ----------
        gene : str
            Gene to run the analysis on.
        celltype_col : str
            The column in the cell_info data frame inside coassay_matrix.h5
            containing the cell groupings.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'a')
        if 'tile_significance/' + celltype_col + '/' + gene in f.keys():
            f.close()
            return
        z_d = self.compute_gene_tile_significance_shap(gene, celltype_col)
        p_d = self.compute_gene_tile_significance_signed_rank(gene, celltype_col, z_d)
        write_significance(f, "tile_significance/" + celltype_col + '/' + gene, z_d, p_d)
        f.close()
        return
    
    def compute_gene_tile_significance_shap(self, gene, celltype_col):
        """For each tile in the tile matrix centered around the gene,
        compute Shapley values for each cell cluster in held-out 
        test set and then standardize the Shapley values.
        
        Parameters
        ----------
        gene : str
            Gene to run the analysis on.
        celltype_col : str
            The column in the cell_info data frame inside coassay_matrix.h5
            containing the cell groupings.

        Returns
        -------
        z_d
            Standardized Shapley values.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        w_mat = np.array(f['genes/' + gene][:])
        e = np.array([f['genes/' + gene].attrs['intercept']])
        train_alpha = float(f['genes/' + gene].attrs['alpha'])
        gex_train, gex_test = self.get_gex_gene(gene)
        tile_gene_mat_train, tile_gene_mat_test = self.gene_tile_matrix_scaled(gene,
                        normalization_factor='ReadsInTSS')
        x = np.array(tile_gene_mat_train.todense())
        z_d = set_gene_tile_significance_bootstrapped(x, np.ravel(gex_train.todense()), w_mat, e, self.cell_info.iloc[self.train_ix], celltype_col, self.cell_info[celltype_col].unique())
        f.close()
        return z_d 

    def compute_gene_tile_significance_signed_rank(self, gene, celltype_col, z_d):
        """For each tile in the tile matrix centered around the gene,
        perform signed rank test between predicted gene expression and 
        predicted expression without the given tile.
        
        Parameters
        ----------
        gene : str
            Gene to run the analysis on.
        celltype_col : str
            The column in the cell_info data frame inside coassay_matrix.h5
            containing the cell groupings.
        z_d : array
            Standardized Shapley values.

        Returns
        -------
        p_d
            p-values for each tile.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        w_mat = np.array(f['genes/' + gene][:])
        e = np.array([f['genes/' + gene].attrs['intercept']])
        train_alpha = float(f['genes/' + gene].attrs['alpha'])
        gex_train, gex_test = self.get_gex_gene(gene)
        tile_gene_mat_train, tile_gene_mat_test = self.gene_tile_matrix_scaled(gene,
                        normalization_factor='ReadsInTSS')
        x = np.array(tile_gene_mat_test.todense())
        p_d = set_gene_tile_significance_signed_rank(x, np.ravel(gex_test.todense()), w_mat, e, self.cell_info.iloc[self.test_ix], celltype_col, self.cell_info[celltype_col].unique(), z_d)
        f.close()
        return p_d 

    def build_model(self, atac_shape, a):
        """Model framework for run regression.
        
        Parameters
        ----------
        atac_shape : int
            Number of tiles in tile matrix.
        a : float
            Regularization parameter for L2-regularization.

        Returns
        -------
        m
            Model.
        """
        # regression framework
        inputs = tf.keras.layers.Input(shape=(atac_shape,), name = 'inputA')
        out = tf.keras.layers.Dense(1, activation = tf.exp, name = 'rate', kernel_regularizer=tf.keras.regularizers.l2(a), kernel_constraint = tf.keras.constraints.NonNeg())(inputs)
        m = tf.keras.models.Model(inputs=inputs, outputs=out)
        return m

    def get_model_weight_average(self, weights):
        """Compute average regression coefficients and bias estimated
        over each fold of cross-validation.
        
        Parameters
        ----------
        weights : matrix
            Learned regression coefficients and bias across all models trained
            during cross validation. 

        Returns
        -------
        new_weights
            Average regression coefficients and bias.
        """
        new_weights = list()
        for weights_list_tuple in zip(*weights): 
            new_weights.append(np.array([np.array(w).mean(axis=0) for w in zip(*weights_list_tuple)]))
        return new_weights

    
    def run_model_cross_validation(self, rna, atac, epochs, verbose, plot_loss):
        """Perform cross-validation of gene-level regression model.

        Parameters
        ----------
        rna : vector
            Gene expression vector.
        atac : matrix
            Tile-matrix centered on gene body.
        epochs : int
            Number of epochs for model training.
        verbose : bool
            Whether to print progress bar.
        plot_loss : bool
            Whether to plot train and test loss.

        Returns
        -------
        s_best_corr, s_params, best_w
            Average Spearman correlation across all folds for chosen regularization 
            parameter, chosen regularization parameter, average regression coefficients
            of all folds for chosen regularization parameter.
        """

        kfold = 5 
        cv = RepeatedKFold(n_splits=kfold, n_repeats = 1, random_state=9)

        best_corr = -1
        s_best_corr = -1
        best_testloss = 1000000
        flag = 0

        # print("Performing " + str(int(kfold)) + "-fold cross validation...")
        logging.info("Performing " + str(int(kfold)) + "-fold cross validation...")
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
                hist = model_custom.fit(x=atac_train, y=rna_train, validation_data=(atac_test, rna_test), epochs=epochs, verbose=None)
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
            weights = self.get_model_weight_average(model_weights) 
            
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
        """Compute Spearman correlation using the trained
        model on held-out test-set.

        Parameters
        ----------
        rna : vector
            Gene expression vector.
        atac : matrix
            Tile-matrix centered on gene body.
        w : vector
            Learned regression coefficients.
        a : float
            Learned bias.

        Returns
        -------
        s_corr
            Spearman correlation.
        """

        model_custom = self.build_model(atac.shape[1], a)
        model_custom.set_weights(w)
        s_corr, s_pval = self.find_correlation_spearman(model_custom, atac, rna)
        return s_corr

    def check_if_calculated(self, gene):
        """Check if regression model is already trained.

        Parameters
        ----------
        gene : str
            Gene to run regression on.

        Returns
        -------
        cond
            Returns True is already computed.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'a')
        cond = gene in f['genes']
        f.close()
        return cond
    
    def train_test_model(self, gene, normalization_factor='ReadsInTSS', max_zero_fraction=0.9, epochs=20, verbose=True, plot_loss=False, force=False):
        """Perform cross validation for each regularization parameter
        and then choose the best model and test model performance on 
        held-out test set.

        Parameters
        ----------
        gene : str
            Gene to run regression on.
        normalization_factor : str
            Normalization factor of tile matrix.
        max_zero_fraction : float
            Maximum allowed sparsity in gene expression vector.
        epochs : float
            Number of epochs of model training.
        verbose : bool 
            Whether to print progress bar.
        plot_loss : bool
            Whether to plot train and test loss.
        force : bool
            Not implemented. Whether to retrain the gene-model.
        """

        # print("Training regression model on " + gene)
        logging.info("Training regression model on " + gene)
        if self.check_if_calculated(gene) and not force:
            # print("Gene regression model for " + gene + " already calculated. Returning...")
            logging.info("Gene regression model for " + gene + " already calculated. Returning...")
            return
        gex_train, gex_test = self.get_gex_gene(gene)
        gex_train_sparsity = self.check_gex_sparsity(gex_train)
        if gex_train_sparsity > max_zero_fraction:
            # print(gene + " expression too sparse with sparsity = " + str(gex_train_sparsity) + ". Maximum number of zeros allowed is " + str(max_zero_fraction) + ". Returning...")
            logging.info(gene + " expression too sparse with sparsity = " + str(gex_train_sparsity) + ". Maximum number of zeros allowed is " + str(max_zero_fraction) + ". Returning...")
            return

        tile_gene_mat_train, tile_gene_mat_test = self.gene_tile_matrix_scaled(gene, normalization_factor)
        train_corr, train_alpha, train_w = self.run_model_cross_validation(np.ravel(gex_train.todense()), tile_gene_mat_train, epochs, verbose, plot_loss)
        # print("Avg. cross-validation spearman correlation", train_corr)
        # print("Chosen regularization parameter:", train_alpha) 
        logging.info("Avg. cross-validation spearman correlation: " + str(train_corr))
        logging.info("Chosen regularization parameter: " + str(train_alpha))
        if train_corr is None:
            # print("ERROR: Regression could not be estimated. Returning...")
            logging.info("ERROR: Regression could not be estimated. Returning...")
            return
        test_corr = self.test_model(np.ravel(gex_test.todense()), tile_gene_mat_test, train_w, train_alpha)
        # print("Spearman corr on test set:", test_corr)
        logging.info("Spearman corr on test set: " + str(test_corr))

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
        dset.attrs['tilesize'] = tiles.iloc[0]['end'] - tiles.iloc[0]['start']
        dset.attrs['max_zero_fraction'] = max_zero_fraction
        dset.attrs['epochs'] = epochs
        dset.attrs['scATAC norm factor'] = normalization_factor
        f.close()


    def get_gene_corr(self, gene):
        """Get Spearman correlation on held-out test set
        for given gene.

        Parameters
        ----------
        gene : str
            Gene to get correlation for.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            if self.log_file_name == '': print("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            else: logging.info("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            f.close()
            return
        d = {}
        d['spearman_correlation_train'] = f['genes/' + gene].attrs['spearman_correlation_train']
        d['spearman_correlation_test'] = f['genes/' + gene].attrs['spearman_correlation_test']
        return d
    
    def get_gene_coefficient(self, gene):
        """Get regression coefficients for gene-level model.

        Parameters
        ----------
        gene : str
            Gene to get regression coefficients for.

        Returns
        -------
        w_mat, intercept
           Regression coefficients and bias.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            if self.log_file_name == '': print("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            else: logging.info("ERROR: " + gene + " not present in " + self.output_dir + self.out_file) 
            f.close()
            return
        w_mat = f['genes/' + gene][:]
        intercept = f['genes/' + gene].attrs['intercept']
        f.close()
        return w_mat, intercept

    def get_gene_alpha(self, gene):
        """Get chosen regularization parameter for given gene.

        Parameters
        ----------
        gene : str
            Gene to get regularization parameter for.

        Returns
        -------
        alpha
           Regularization parameter.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            if self.log_file_name == '': print("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            else: logging.info("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            f.close()
            return
        alpha = f['genes/' + gene].attrs['alpha']
        f.close()
        return alpha

    def get_gene_window_coords(self, gene):
        """Get coordinates of gene window.

        Parameters
        ----------
        gene : str
            Gene to get coordinates for.

        Returns
        -------
        chrm, start, end
           Gene window coordinates.
        """

        f = h5py.File(self.output_dir + self.out_file, mode = 'r')
        if gene not in f['genes'].keys():
            if self.log_file_name == '': print("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            else: logging.info("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            f.close()
            return
        chrm = f['genes/' + gene].attrs['chr']
        start = int(f['genes/' + gene].attrs['start'])
        end = int(f['genes/' + gene].attrs['end'])
        return chrm, start, end
    
    def plot_gene(self, gene, groups = 'Clusters', plot_frags = False, to_save = False, plot_dir='', cmap = None, save_format='png', figsize=(17, 14), sort_gex=False, show_yticks=False, plot_shap=False, shap_cmap='Blues', cluster_order=[], bg_transparent=False):

        """Plot SCARlink output for given gene.

        Parameters
        ----------
        gene : str
            Gene to plot output for.
        groups : str
            Cell grouping label. This should be present as a column name
            in the data frame cell_info found in coassay_matrix.h5.
        plot_frags : bool
            Not fully implemented. Whether to plot accessibility at base-pair
            resolution or at tile-level.
        to_save : bool
            Whether to save output plot.
        plot_dir : str
            Path to directory where output plot is to be saved. plot_dir='' 
            will save plot in present working directory.
        cmap : dictionary
            Color map for each cell cluster.
        save_format : str
            File format for saving output plot.
        figsize : (float, float)
            Figure size for output plot.
        sort_gex : bool
            Whether to order the cell clusters based on average expression
            of gene.
        show_yticks : bool
            Whether to show yticks for accessibility in tiles.
        plot_shap : bool
            Whether to plot standardized Shapley values for each cell cluster.
        shap_cmap : str or colormap
            Color map for plotting Shapley values.
        cluster_order : list
            Order in which to plot clusters.
        bg_transparent : bool
            Whether to make background of output plot transparent.
        """

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
            if self.log_file_name == '': print("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            else: logging.info("ERROR: " + gene + " not present in " + self.output_dir + self.out_file)
            f.close()
            return
        w_mat = f['genes/' + gene][:] 
        chrm = f['genes/' + gene].attrs['chr']
        start = int(f['genes/' + gene].attrs['start'])
        end = int(f['genes/' + gene].attrs['end'])
        tilesize = int(f['genes/' + gene].attrs['tilesize'])
        sp_corr = f['genes/' + gene].attrs['spearman_correlation_test']

        if plot_shap:
            zscore_d = self.compute_gene_tile_significance_shap(gene, groups)
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

        if plot_shap: 
            ax_shap = plt.subplot(gs[0, 1])
            s_plot = ax_shap.scatter([0, 0], [0, 0], cmap=shap_cmap, vmin=min_pvals, vmax=max_pvals, c=[0, 0])
            ax_shap.set_xlim((1, 2))
            ax_shap.set_ylim((1, 2))
            ax_shap.axis('off')
            cb = plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min_pvals, vmax=max_pvals), cmap=shap_cmap), ax=ax_shap, orientation='horizontal', label='z-score', location='top')
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
            
            if plot_shap:
                max_pvals = 4
                pvals_d = zscore_d
                ylm = ax[c][0].get_ylim()
                yrange = ylm[1] - ylm[0]
                # ordering so that higher absolute shap values are plotted last
                o_idx = np.argsort(zscore_d[clust])
                ax[c][0].scatter(x=plotrange[o_idx], y=[ylm[0] - 0.1*yrange]*zscore_d[clust].shape[0], marker='o', alpha=0.5, clip_on=False, c=pvals_d[clust][o_idx], cmap=shap_cmap, vmin=0, vmax=max_pvals)
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

            if not show_yticks: 
                ax[c][0].set_yticks([])
                ax[c][0].spines['left'].set_visible(False) 
            
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

        # plt.suptitle(gene + "(corr : " + str(round(sp_corr, 4)) + ")", fontsize = 'xx-large')
        plt.suptitle(gene + "(" + r'$\rho$: ' + str(round(sp_corr, 4)) + ")", fontsize = 'xx-large')
        if to_save:
            plt.savefig(plot_dir + gene + '.' + save_format, transparent=bg_transparent)
            if self.log_file_name == '': print("Saved as " + plot_dir + gene + '.' + save_format) 
            else: logging.info("Saved as " + plot_dir + gene + '.' + save_format) 
        
