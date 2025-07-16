#!/bin/bash
#SBATCH --job-name=scRNA-analysis
#SBATCH --partition=msc_appbio
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/scRNA_prep/scripts/slurms/scRNA_analysis/scRNA_alignment-%j.out

# Directory for input files

BASE_DIR="/scratch/prj/bmb_tofacitinib"
BASE_DATA_DIR="${BASE_DIR}/data"

BASE_SINGULARITY_DIR="${BASE_DIR}/research_project/singularity_images"
SINGULARITY_IMG_DIR="${BASE_SINGULARITY_DIR}/cdseqv1.0.9.sif"

OUTPUT_DIR="${BASE_DATA_DIR}/outputs/scRNA_processed"
mkdir -p $OUTPUT_DIR


singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind ${BASE_DIR}:${BASE_DIR} \
  $SINGULARITY_IMG_DIR Rscript cibersortx_prep.R

echo "CibersortX Single Cell Reference Formatting Complete (CDseq v1.0.9)."
