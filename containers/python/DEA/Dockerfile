
FROM python:3.13.2-bookworm

RUN apt-get update && apt-get install -y wget bzip2 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh



# Adding conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Initialize Conda for bash
RUN conda init bash

# Create the pydeseq2 environment based on version python 3.12.2 
RUN conda create -n pydeseq2_env -c conda-forge -c bioconda pydeseq2 python=3.13.2


# Automatically activate the conda environment when running iteratively to test libraries
RUN echo "conda activate pydeseq2_env" >> ~/.bashrc
