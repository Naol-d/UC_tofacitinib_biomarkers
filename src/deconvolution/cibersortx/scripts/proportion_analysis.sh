#!/bin/bash
#SBATCH --job-name=cell_proportion_analysis
#SBATCH --partition=msc_appbio
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=7G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/deconvolution/cibersortx/scripts/slurms/proportion_analysis/cell_prop_analysis-%j.out
#SBATCH --time=48:00:00

# Base Directories
BASE_DIR="/scratch/prj/bmb_tofacitinib"
SINGULARITY_IMG_DIR="${GENERAL_ANALYSIS_DIR}/singularity_images/cdseqv1.0.9.sif"

singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind $BASE_DIR:$BASE_DIR \
  $SINGULARITY_IMG_DIR Rscript cell_proportion_analysis.R

echo "Cibersortx Fractions Analysis Complete."
