from setuptools import setup
import platform

config = {
    'name': 'scarlink',
    'description': 'Single-cell ATAC+RNA linking using regularized Poisson regression',
    'version': '1.0.0',
    'packages': ['scarlink'],
    'setup_requires': [],
    'install_requires': ['shap==0.41.0',
                         'tables==3.7.0',
                         'scikit-learn==1.2.1',
                         'numpy==1.23.1',
                         'scanpy==1.9.3',
                         'fa2',
                         'python-igraph==0.9.11',
                         'rpy2>=3.5.11'] + ([] if platform.platform().startswith('macOS') else ['tensorflow==2.11.0']),
    'entry_points': {'console_scripts': [
        'scarlink = scarlink.run_scarlink:main',
        'scarlink_tiles = scarlink.run_scarlink_tiles:main',
        'scarlink_plot = scarlink.run_scarlink_visualization:main',
        'scarlink_processing = scarlink.preprocessing.create_h5_files:main',
        'scarlink_add_cell_cluster = scarlink.run_scarlink_add_cluster:main']}
}

if __name__== '__main__':
    setup(**config)
