#!/bin/bash
#SBATCH --job-name=deconvolution_cibersortxfractions
#SBATCH --partition=msc_appbio
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/deconvolution/cibersortx/scripts/slurms/cibersortxfractions/fractions-%j.out
#SBATCH --time=48:00:00

# Base Directory
GENERAL_ANALYSIS_DIR="/scratch/prj/bmb_tofacitinib/research_project"
BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"
BASE_DIR="${GENERAL_ANALYSIS_DIR}/deconvolution/cibersortx"

# Directory for input and output files
INPUT_DIR="${BASE__DATA_DIR}/outputs/bulkRNA_deconvolution_prep"
OUTPUT_DIR="${BASE__DATA_DIR}/outputs/deconvolution/cibersortx/TPM_fractions_0.5"
mkdir -p $OUTPUT_DIR

echo "Currently IP Address..."

curl https://ifconfig.me

echo ""

SINGULARITY_IMG_DIR="${GENERAL_ANALYSIS_DIR}/singularity_images/fractions_latest.sif"
CIBERSORTX_TOKEN="${GENERAL_ANALYSIS_DIR}/singularity_images/tokens/cibersortx_token.txt"
TOKEN=$(head -n 1 "$CIBERSORTX_TOKEN")

singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind $INPUT_DIR:/src/data \
  --bind $OUTPUT_DIR:/src/outdir \
  $SINGULARITY_IMG_DIR /src/CIBERSORTxFractions --username k24103029@kcl.ac.uk --token $TOKEN \
	--single_cell TRUE \
	--mixture /src/data/corrected_TPM_counts_bulkRNA_CSX.txt \
	--refsample /src/data/integrated_RKPM_count_CIBERSORTx_reference.txt \
	--rmbatchSmode TRUE \
	--perm 1000 \
	--fraction 0.5 \
	--outdir /src/outdir \
	--verbose TRUE

echo "Cibersortx Fractions Analysis Complete (Cibersortx Fractions)."
