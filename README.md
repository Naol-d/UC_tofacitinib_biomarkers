# Identifying markers of Relapse in Ulcerative Colitis Patients when Tofacitinib dosage is reduced
The pipeline starts from alignment to statisitical analysis uncovering differentially expressed genes and changes in cell populations across bulkRNA-seq data through deconvolution techniques.

## Reproducibility
### Utilisation of Docker Images in a HPC or Multi-User Platform
Docker images are available on the DockerHub Repository for direct pulling. Docker images can then be converted to Singularity Image Formats (.sif).  
**Example Command:** singularity pull docker://<"name of docker image">

**Access Docker images using commands:**
```bash
docker pull naolxi/cellranger:9.0.1
docker pull naolxi/pydeseq2:0.5.0
docker pull naolxi/r-env:4.4.1
docker pull naolxi/deeplearn:1
```
**Conda environments** utilised can be imported using .yml files available in the 'conda_envs' folder 

## How to Run
Pipeline is modular, bash scripts have been provided for HPC friendly environments based on SLURM job scheduling. The order of execution can be identified in the Directory Structure 'src' subtree section below.

## Directory Structure
```text
UC_tofacitinib_biomarkers/
├── .gitignore
├── LICENSE
├── README.md
├── conda_envs/
│   ├── .gitkeep
│   ├── sratools_env.yml                        # Subread (featureCounts) environment for portability
│   └── subread_env.yml                         # SRA tools environment for portability
├── containers/
│   ├── .gitkeep
│   ├── cell_ranger/                            # Dockerfile for cell ranger container
│   ├── python/                                 # Dockerfile for Python‑based DEA & ML
│   └── R/                                      # Dockerfile for R environment used in RNA analysis
├── data/                                       # Data not publicly available (access on publication)
│   └── .gitkeep
└── src/
    ├── bulkRNA_alignment/                      # Bulk RNA alignment
    │   └── scripts/
    │       ├── batch1_disease_alignment.sh
    │       ├── batch2_disease_alignment.sh
    │       ├── featureCounts.sh
    │       ├── generate_genome_indices.sh
    │       └── slurms/
    ├── bulkRNA_prep/                           # Matrix prep, batch correction, GSEA formatting
    │   └── scripts/
    │       ├── bulkRNA_prep.sh
    │       └── bulkRNA_prep.R
    ├── bulkRNA_DEA/                            # DEA, GO & KEGG enrichment
    │   └── scripts/
    │       ├── analysis.sh
    │       └── DEG_analysis.R
    ├── scRNA_alignment/                        # Retrieve & align scRNA data
    │   ├── download/
    │   │   ├── dataset_1_ssr_list.txt
    │   │   ├── dataset1_download.sh
    │   │   ├── dataset_2_srr_list.txt
    │   │   └── dataset_2_download.sh
    │   └── alignment/
    │       ├── build_ref.sh
    │       ├── rename_dataset_1_fastqs.sh
    │       ├── rename_dataset_2_fastqs.sh
    │       ├── scRNA_dataset1_alignment.sh
    │       └── scRNA_dataset2_alignment.sh
    ├── scRNA_prep/                             # Prep scRNA reference matrix for signatures
    │   └── scripts/
    │       ├── cibersortx_prep.R
    │       └── cibersortx_prep.sh
    └── deconvolution/                          
        └── cibersortx/                        # bulkRNA deconvolution & cell‑type analysis
            └── scripts/
                ├── browser/                   # (Optional) Allows for web‑browser capabilities in a HPC environment required to generate a cibersortX token linked to a compute node's IP
                │   └── browser.sh
                ├── deconvolution_bulkRNA_prep.R
                ├── deconvolution_bulkRNA_prep.sh
                ├── cibersortx_fractions.sh
                ├── cell_proportion_analysis.R
                └── proportion_analysis.sh

``` 



