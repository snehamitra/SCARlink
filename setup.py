from setuptools import setup

config = {
    'name': 'scarlink',
    'description': 'Single-cell ATAC+RNA linking using regularized Poisson regression',
    'version': '0.0.1',
    'packages': ['scarlink'],
    'setup_requires': [],
    'install_requires': ['tensorflow==2.11.0',
                         'shap==0.41.0',
                         'tables==3.7.0',
                         'scikit-learn==1.2.1',
                         'numpy==1.23.1',
                         'scanpy==1.9.3',
                         'fa2',
                         'python-igraph==0.9.11',
                         'rpy2==3.5.12'],
    'entry_points': {'console_scripts': [
        'scarlink = scarlink.run_scarlink:main',
        'scarlink_tiles = scarlink.run_scarlink_tiles:main',
        'scarlink_plot = scarlink.run_scarlink_visualization:main',
        'scarlink_processing = scarlink.preprocessing.create_h5_files:main']}
}

if __name__== '__main__':
    setup(**config)
