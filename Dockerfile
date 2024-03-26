FROM ubuntu:22.04

# Set the working directory inside the container
WORKDIR /app

# Install basic utilities
RUN apt-get -y update
RUN apt-get -y install bzip2 ca-certificates curl unzip wget

# Download and install miniconda
RUN curl -L https://repo.anaconda.com/miniconda/Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -O && bash Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -bf -p /opt/conda && rm Miniconda3-py39_23.11.0-2-Linux-x86_64.sh

# Set conda path
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && echo "conda activate base" >> ~/.bashrc 
ENV PATH /opt/conda/bin:$PATH

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
RUN pip install notebook
