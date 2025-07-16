#!/bin/bash
#SBATCH --job-name=human_genome_index
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_alignment/scripts/slurms/diseased_alignments/generate_genome/slurm-%j.out

# Load the STAR module
module load star/2.7.10b-gcc-13.2.0

BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"

REFERENCE_GENOME_DIR="${BASE_DATA_DIR}/reference_GRCh38_gencode_v47"

BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_alignment/human_genome"

GENOME_FASTA="${REFERENCE_GENOME_DIR}/GRCh38.primary_assembly.genome.fa"
GTF_FILE="${REFERENCE_GENOME_DIR}/gencode.v47.primary_assembly.annotation.gtf"

GENOME_OUTPUT="${BASE_DATA_DIR}/outputs/bulkRNA_aligned_outputs/human_genome_index"
mkdir -p $GENOME_OUTPUT

# Generate the genome indices
STAR --runMode genomeGenerate \
     --genomeDir $GENOME_OUTPUT \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 149 \
     --runThreadN 16
