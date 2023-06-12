import h5py
from scarlink.src.model import RegressionModel

def read_model(out_dir, out_file_name = 'coefficients.hd5', input_file_name='', read_only=False):
    # read_only=True: mode = 'r', else mode = 'a'
    out_dir = out_dir + '/' if out_dir[-1] != '/' else out_dir
    
    f = h5py.File(out_dir + out_file_name, mode = 'r')
    # if input_file_name == '': input_file_name = f['genes'].attrs['input_file_name']
    dirname = '/'.join(out_dir[:-1].split('/')[:-1]) + '/'
    input_file_name = dirname + 'coassay_matrix.h5'
    gtf_file = f['genes'].attrs['gtf_file']
    scatac_fragment_file = f['genes'].attrs['scatac_fragment_file']
    f.close()

    if read_only: mode = 'r'
    else: mode = 'a'
    m = RegressionModel(input_file = input_file_name, output_dir = out_dir, gtf_file = gtf_file, scatac_fragment_file = scatac_fragment_file, out_file_name = out_file_name, mode=mode)
    return m
