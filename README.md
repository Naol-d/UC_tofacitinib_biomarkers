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
**Conda environments** utilised used be imported using .yml files available in the 'conda_envs' folder 

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
│   ├── python/                                 # Dockerfile to faciliate python alternative DEA; Machine Learning
│   └── R/                                      # Dockerfile for R Environment used in RNA analysis
├── data/                                       # Data is currently not publicly available (Access will available upon publication)
│   └── .gitkeep
└── src/
    ├── bulkRNA_alignment/                      # BulkRNA alignment 
    │   └── scripts/
    │       ├── batch1_disease_alignment.sh
    │       ├── batch2_disease_alignment.sh
    │       ├── featureCounts.sh
    │       ├── generate_genome_indices.sh
    │       └── slurms/
    ├── bulkRNA_prep/                           # Preprocessing of matrix data, batch correction, GSEA input pre-formatting
    │   └── scripts/
    │       ├── bulkRNA_prep.sh
    │       └── bulkRNA_prep.R
    ├── bulkRNA_DEA/                            # Differential Expression Analysis, GO Enrichment, KEGG Pathway Enrichment
    │   └── scripts/
    │       ├── analysis.sh
    │       └── DEG_analysis.R
    ├── scRNA_alignment/                        
    │   ├── download/                           # Retrieval of public scRNA datasets
    │   │   ├── dataset_1_ssr_list.txt
    │   │   ├── dataset1_download.sh
    │   │   ├── dataset_2_srr_list.txt
    │   │   └── dataset_2_download.sh
    │   └── alignment/                          # File namming format processing and scRNA alignment
    │       ├── build_ref.sh
    │       ├── rename_dataset_1_fastqs.sh
    │       ├── rename_dataset_2_fastqs.sh
    │       ├── scRNA_dataset1_alignment.sh
    │       └── scRNA_dataset2_alignment.sh
    ├── scRNA_prep/                             # Preparation of scRNA reference matrix file for signature matrix generation (single cell processing)
    │   └── scripts/
    │       ├── cibersortx_prep.R
    │       └── cibersortx_prep.sh
    └── deconvolution/
        └── cibersortx/                         # bulkRNA deconvolution using CIBERSORTX and cell population analysis  
            ├── cell_proportion_analysis.R
            ├── cibersortx_fractions.sh
            ├── corrected_bulkRNA_prep.R
            ├── proportion_analysis.sh
            └── slurms/
``` 



