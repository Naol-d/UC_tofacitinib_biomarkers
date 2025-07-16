#!/bin/bash
#SBATCH --job-name=bulkRNA_prepping
#SBATCH --cpus-per-task=1
#SBATCH --partition=msc_appbio
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_prep/scripts/slurms/bulk_prep-%j.out

# Directory for input files

BASE_DIR="/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_prep"
BASE_SINGULARITY_DIR="/scratch/prj/bmb_tofacitinib/research_project/singularity_images"

INPUT_DIR="${BASE_DIR}/originals"

R_SINGULARITY_IMG_DIR="${BASE_SINGULARITY_DIR}/cdseqv1.0.9.sif"

singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind ${BASE_DIR}:${BASE_DIR} \
  $R_SINGULARITY_IMG_DIR Rscript bulkRNA_prep.R

echo "Preliminary Data Analysis and Counts Data Batch Correction Complete"

