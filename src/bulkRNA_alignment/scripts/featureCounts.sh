#!/bin/bash
#SBATCH --job-name=gene_quantification
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_alignment/scripts/slurms/diseased_alignments/feature_counts/slurm-%j.out

BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"

# Directory for aligned reads
ALIGNED_DIR="${BASE_DATA_DIR}/outputs/bulkRNA_aligned_outputs"

# Directory containing the human genome annotations
HUMAN_GENOME_ANNON_DIR="${BASE_DATA_DIR}/originals/reference_GRCh38_gencode_v47/gencode.v47.primary_assembly.annotation.gtf"

eval "$(conda shell.bash hook)"
conda activate subread

OUTPUT_DIR="/scratch/prj/bmb_tofacitinib/data/outputs/bulkRNA_aligned_outputs"
mkdir -p $OUTPUT_DIR

featureCounts -T 8 -p --countReadPairs -t exon -g gene_id -a $HUMAN_GENOME_ANNON_DIR -o ${OUTPUT_DIR}/FC_combo.txt ${ALIGNED_DIR}/batch_2_disease/*.bam ${ALIGNED_DIR}/batch_1_disease/*.bam
