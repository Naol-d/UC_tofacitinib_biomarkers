#!/bin/bash
#SBATCH --job-name=DEG_analysis
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/DEG/scripts/slurms/analysis-%j.out

# Directory for input files

BASE_DIR="/scratch/prj/bmb_tofacitinib"
BASE_SINGULARITY_DIR="/scratch/prj/bmb_tofacitinib/research_project/singularity_images"

INPUT_DIR="${BASE_DIR}/originals"


R_SINGULARITY_IMG_DIR="${BASE_SINGULARITY_DIR}/cdseqv1.0.9.sif"


singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind ${BASE_DIR}:${BASE_DIR} \
  $R_SINGULARITY_IMG_DIR Rscript DEG_analysis.R

echo " DEG Analysis Done..."
