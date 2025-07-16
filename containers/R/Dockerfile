# Stage 1:
FROM kitware/cmake:ci-debian12-x86_64-2025-01-23 AS cmake-build

RUN apt-get update && apt-get install -y \
    libgsl-dev


# Stage 2
# R 4.4.1
FROM rocker/r-ver:4.4.1

# Copy compiled dependencies from the CMake stage, libgsl-dev For RcppGSL gets copied to the second stage build
COPY --from=cmake-build /usr/lib /usr/lib
COPY --from=cmake-build /usr/bin /usr/bin
COPY --from=cmake-build /usr/include /usr/include

# Installing system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    ghostscript \
    libglu1-mesa-dev \
    libgl1-mesa-dev \
    mesa-common-dev \
    libosmesa6-dev \
    freeglut3-dev \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    libhiredis-dev \
    && rm -rf /var/lib/apt/lists/*


# Add OpenMPI library path to LD_LIBRARY_PATH so dynamic linker can find MPI shared libraries
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib:$LD_LIBRARY_PATH

# Install BiocManager 
RUN R -e 'install.packages("BiocManager", repos="http://cran.r-project.org", dependencies=TRUE)'

# Initalise BiocManager
RUN R -e 'BiocManager::install()'


# No or Base R Dependencies Requiring Packages
RUN R -e 'install.packages(c( \
    "Rcpp", "MASS", "foreach", "doParallel", \
    "dirmult", "RcppThread", "iterators", \
    "magrittr", "rlang", "matrixStats", \
    "devtools", "compositions", "doRNG", "randtests" \   
    ), repos="http://cran.r-project.org", dependencies=TRUE, quiet=TRUE)'


RUN R -e 'BiocManager::install(c( \
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", \
    "SummarizedExperiment", "SingleCellExperiment", "DESeq2", "MAST", \
    "monocle", "limma", "DelayedArray", "multtest", "rtracklayer", "sva", \
    "variancePartition", "org.Hs.eg.db"))'


# Dependencies Requiring Packages
RUN R -e 'install.packages(c( \
    "Seurat", "ggplot2", "dplyr", "Matrix", \
    "ggpubr", "harmony", "pheatmap" \
    ),repos="http://cran.r-project.org", dependencies=TRUE)'


# Install CDSeq from GitHub
RUN R -e 'devtools::install_github("kkang7/CDSeq_R_Package")'

# Installing ggtangle needed for enrichplot used in clusterProfiler and harmony for batch correction
RUN R -e 'install.packages(c("ggtangle", "harmony"), repos="http://cran.r-project.org", dependencies=TRUE)'

# Install clusterProfiler
# Install biomaRt, goseq and scDblFinder
RUN R -e 'BiocManager::install(c("clusterProfiler","goseq", "biomaRt", "scDblFinder", "glmGamPoi"), ask = FALSE)'

# Install presto from GitHub (for more efficient compute in marker identification)
RUN R -e 'devtools::install_github("immunogenomics/presto")'

RUN R -e 'install.packages(c("car", "vegan", "effsize", "DirichletReg"), repos="http://cran.r-project.org", dependencies=TRUE)'

# Assessing batch effects in scRNA data
RUN R -e 'devtools::install_github("theislab/kBET")'


