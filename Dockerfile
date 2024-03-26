FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Install basic utilities
RUN apt-get -y update
RUN apt-get -y install bzip2 ca-certificates 

# Create conda env
RUN conda create -n scarlink-env python=3.8

# Override default shell and use bash
SHELL ["conda", "run", "-n", "scarlink-env", "/bin/bash", "-c"]

# Set conda channel order
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install packages from conda
RUN conda install -y -c conda-forge mamba
RUN mamba install -y -c conda-forge r-seurat r-devtools r-biocmanager 
RUN mamba install -y -c bioconda bioconductor-rhdf5 bioconductor-genomeinfodbdata bioconductor-chromvar bioconductor-motifmatchr bioconductor-complexheatmap
RUN mamba install -y -c conda-forge 'rpy2>=3.5.11'
RUN mamba install -y -c conda-forge fa2

# Install ArchR
RUN R -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())'

# Set environment variables for Python
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Install SCARlink
RUN mkdir /app/SCARlink
COPY . /app/SCARlink
WORKDIR /app
RUN pip install -e SCARlink

