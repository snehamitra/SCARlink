#################################################################################
# Read contents of h5 file and group cells if group_cells = True
#################################################################################

import h5py
import pandas
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.model_selection import train_test_split

#def add_df_rows(rows):
    
def construct_cell_info(f, group_cells):
    cell_info = f.select('cell_info') # _annotated')
    if not group_cells:
        return cell_info
    train_cells = f.get_node('group_cells')['train'][:]
    test_cells = f.get_node('group_cells')['test'][:]
    cells = np.append(train_cells, test_cells, axis = 1)

    group_cell_info = pandas.DataFrame(columns = list(cell_info))
    for col in list(cell_info):
        if cell_info[col].dtype != 'object':
            group_cell_info[col] = cell_info[col].values[np.ravel(cells)].reshape(cells.shape).sum(axis = 0)
            group_cell_info = group_cell_info.astype({col: cell_info[col].dtype})
                
        else:
            group_cell_info[col] = pandas.DataFrame(cell_info[col].values[np.ravel(cells)].reshape(cells.shape)).T.mode(axis = 1)[0]

    # update cell names
    cell_name = ['train_' + str(i) for i in range(train_cells.shape[1])] + ['test_' + str(i) for i in range(test_cells.shape[1])]
    group_cell_info['cell_name'] = cell_name
    return group_cell_info

def construct_gex_mat(f, cell_info, group_cells):
    gex_matrix = csr_matrix((f.get_node('gene_expression')['data'][:], f.get_node('gene_expression')['indices'][:], f.get_node('gene_expression')['indptr'][:]), shape=np.flip(f.get_node('gene_expression')['shape'][:]))
    if not group_cells:
        return gex_matrix
    gex_matrix = csr_matrix((f.get_node('gene_expression_raw')['data'][:], f.get_node('gene_expression_raw')['indices'][:], f.get_node('gene_expression_raw')['indptr'][:]), shape=np.flip(f.get_node('gene_expression_raw')['shape'][:]))

    train_cells = f.get_node('group_cells')['train'][:]
    test_cells = f.get_node('group_cells')['test'][:]
    cells = np.append(train_cells, test_cells, axis = 1)
    group_gex_matrix = gex_matrix[cells[0]]
    for i in range(1, cells.shape[0]):
        group_gex_matrix += gex_matrix[cells[i]]

    # scale values
    row_indices, col_indices = group_gex_matrix.nonzero()
    scale_factor = (1e6 / cell_info['nCount_RNA'].values)[row_indices]
    group_gex_matrix.data *= scale_factor
    return group_gex_matrix

def get_train_test_split(f, gex_matrix, random_state, group_cells):
    train_ix, test_ix = train_test_split(np.arange(gex_matrix), random_state = 9)
    if not group_cells:
        return train_ix, test_ix
    train_cells = f.get_node('group_cells')['train'][:]
    test_cells = f.get_node('group_cells')['test'][:]
    group_train_ix = np.arange(train_cells.shape[1])
    group_test_ix = train_cells.shape[1] + np.arange(test_cells.shape[1])
    return group_train_ix, group_test_ix

def get_gene_tile_matrix_group_cells(f, gene, group_cells):
    tile_gene_mat = csr_matrix((f.get_node(gene)['data'][:], f.get_node(gene)['indices'][:], f.get_node(gene)['indptr'][:]), shape=np.flip(f.get_node(gene)['shape'][:]))
    if not group_cells: return tile_gene_mat

    train_cells = f.get_node('group_cells')['train'][:]
    test_cells = f.get_node('group_cells')['test'][:]
    cells = np.append(train_cells, test_cells, axis = 1)
    group_tile_gene_mat = tile_gene_mat[cells[0]]
    for i in range(1, cells.shape[0]):
        group_tile_gene_mat += tile_gene_mat[cells[i]]
    return group_tile_gene_mat

# take dictionary of z-scores/p-values and create sparse data frame
def sparsify_df(d, logneg=False):
    clusters = sorted(list(d.keys()))
    df = pandas.DataFrame(columns=clusters)
    for ky in d:
        df[ky] = d[ky]
    if logneg:
        df = np.log10(df).abs()
    else:
        df[df < 0] = 0

    return df.columns, csr_matrix(df.values)

def write_sparse_significance(f, df, k):
    g = f.create_group(k)
    g.create_dataset('data', data=df.data)
    g.create_dataset('indptr', data=df.indptr)
    g.create_dataset('indices', data=df.indices)
    g.attrs['shape'] = df.shape

def read_sparse_significance(f, k, entry):
    g = f[k + '/' + entry]
    m = csr_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), g.attrs['shape'])
    return m

# write z-scores and p-values in sparse format
def write_significance(f, k, z, p):
    del f[k]
    g = f.create_group(k)
    clusters, df_z = sparsify_df(z)
    clusters, df_p = sparsify_df(p, logneg=True)

    write_sparse_significance(f, df_z, k+'/z-score')
    write_sparse_significance(f, df_p, k+'/p-value')
    
    sparse_z = read_sparse_significance(f, k, entry='z-score')
    sparse_p = read_sparse_significance(f, k, entry='p-value')

