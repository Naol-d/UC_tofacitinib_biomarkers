#!/bin/bash
#SBATCH --job-name=diseased_alignments
#SBATCH --array=1-10
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/bulkRNA_alignment/scripts/slurms/diseased_alignments/batch_1/slurm-%A_%a.out


echo $SLURM_ARRAY_TASK_ID
# Gets the line position of a fastq name within the disease_samples_batch_1.txt based on the array iteration
Filename=`head disease_samples_batch_1.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1`
echo $Filename

BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"
HUMAN_GENOME_DIR="${BASE_DATA_DIR}/outputs/bulkRNA_aligned_outputs/human_genome_index"

# Directory for all of the batch_1 disease samples
INPUT_DIR="${BASE_DATA_DIR}/originals/bulkRNA/bulkRNA_raw_data/disease/batch_1/"

R1="${INPUT_DIR}/${Filename}_1.fq.gz"
R2="${INPUT_DIR}/${Filename}_2.fq.gz"

# Directory for outputs
OUTPUT_DIR="${BASE_DATA_DIR}/outputs/bulkRNA_aligned_outputs/batch_1_disease/"
mkdir -p $OUTPUT_DIR


module load star/2.7.10b-gcc-13.2.0

STAR --genomeDir $HUMAN_GENOME_DIR \
--runThreadN 16 \
--readFilesIn $R1 $R2 \
--readFilesCommand zcat \
--outFileNamePrefix ${OUTPUT_DIR}/${Filename}_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
