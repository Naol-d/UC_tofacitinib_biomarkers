#!/bin/bash
#SBATCH --job-name=ref_index
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/scRNA_alignment/scripts/slurms/build_ref_genome/build_ref-%j.out

# Directory for input files
BASE_DIR="/scratch/prj/bmb_tofacitinib"
BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"
BASE_SINGULARITY_DIR="/scratch/prj/bmb_tofacitinib/research_project/singularity_images"

CELL_RANGER_DIR="${BASE_SINGULARITY_DIR}/cellrangerv9.0.1.sif"

GENOME_DIR="${BASE_DATA_DIR}/originals/reference_GRCh38_gencode_v47"

OUTPUT_DIR="${BASE_DATA_DIR}/originals/scRNA/scRNA_ref_index"
mkdir -p $OUTPUT_DIR

singularity exec --env LANG=C.UTF-8 --env LC_ALL=C.UTF-8 \
  --bind ${BASE_DIR}:${BASE_DIR} \
  $CELL_RANGER_DIR cellranger mkref \
  --genome=GRCh38_gencode_v47 \
  --fasta=${GENOME_DIR}/GRCh38.primary_assembly.genome.fa \
  --genes=${GENOME_DIR}/gencode.v47.primary_assembly.annotation.gtf \
  --output-dir=$OUTPUT_DIR \
  --nthreads=8

echo "Reference Built (Cell Ranger v9.0.1)."


