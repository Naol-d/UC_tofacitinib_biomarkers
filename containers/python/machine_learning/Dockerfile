FROM python:3.12.10-bookworm

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

# Create the deeplearn environment based on version python 3.12.10 
RUN conda create -n pytorch_env -c conda-forge -c pytorch -c nvidia \
    python=3.12.10 \
    pytorch torchvision torchaudio pytorch-cuda=12.1 \
    jupyterlab notebook ipykernel \
    ipywidgets nb_conda_kernels \
    pandas matplotlib scikit-learn seaborn \
    && conda clean -afy

RUN conda create -n tensorflow_env -c conda-forge \
    python=3.12.10 \
    tensorflow \
    jupyterlab notebook ipykernel \
    ipywidgets nb_conda_kernels \
    pandas matplotlib seaborn scikit-learn \
    && conda clean -afy


# Register the kernel for Jupyter
RUN /bin/bash -c "source activate pytorch_env && \
    python -m ipykernel install --user --name pytorch_env --display-name 'Python (PyTorch)'"

RUN /bin/bash -c "source activate tensorflow_env && \
    python -m ipykernel install --user --name tensorflow_env --display-name 'Python (TensorFlow)'"


# Expose Jupyter port
EXPOSE 8888

# Set a working directory
WORKDIR /workspace

# Use bash as entrypoint to allow dynamic environment switching.
# The default CMD activates pytorch_env and launches Jupyter,
# but this can be overridden at runtime
ENTRYPOINT ["/bin/bash", "-c"]

# Start Jupyter and pytorch_env on container launch
CMD ["/bin/bash", "-c", "source activate pytorch_env && jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root"]
