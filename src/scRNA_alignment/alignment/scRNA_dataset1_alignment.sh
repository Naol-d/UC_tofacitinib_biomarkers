#!/bin/bash
#SBATCH --job-name=scRNA-alignment
#SBATCH --array=1-4
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=70G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/scRNA_alignment/scripts/slurms/scRNA_alignment/dataset1/scRNA_alignment-%A-%a.out

# Directory for input files and container

BASE_DIR="${BASE_DIR}/scratch/prj/bmb_tofacitinib"
BASE_DATA_DIR="${BASE_DIR}/data"
BASE_SINGULARITY_DIR="${BASE_DIR}/research_project/singularity_images"
CELL_RANGER_DIR="${BASE_SINGULARITY_DIR}/cellrangerv9.0.1.sif"


# Where the renamed fastq for each sample exist
FASTQ_DIR="${BASE_DATA_DIR}/originals/scRNA/dataset_1/renamed_fastq"
GENOME_DIR="${BASE_DIR}/originals/scRNA/scRNA_ref_index"
METADATA="${BASE_DIR}/originals/scRNA/metadata/public_metadata.tsv"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
GSM_ID=$(echo "$LINE" | cut -f1)

OUTPUT_DIR="${BASE_DATA_DIR}/outputs/scRNA_aligned_outputs/dataset_2/${GSM_ID}"
mkdir -p $OUTPUT_DIR


singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind ${BASE_DIR}:${BASE_DIR} \
  $CELL_RANGER_DIR cellranger count \
  --id="${GSM_ID}" \
  --create-bam true \
  --transcriptome="${GENOME_DIR}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${GSM_ID}" \
  --localcores=8 \
  --localmem=64 \
  --output-dir="${OUTPUT_DIR}"

echo "Sample ${GSM_ID} Aligned (Cell Ranger v9.0.1)."




