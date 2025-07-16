#!/bin/bash
#SBATCH --job-name=download_dataset1
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/scRNA_alignment/scripts/slurms/data_download/dataset2/public_data_download-%A_%a.out


echo $SLURM_ARRAY_TASK_ID
# Gets the line position of a fastq name within the filename.txt based on the array iteration
SRR=`head dataset_2_srr_list.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1`
echo $SRR

# Activate sra tools
eval "$(conda shell.bash hook)"
conda activate sratools_env

# Where the public scRNA data will be stored
DATA_DIR="/scratch/prj/bmb_tofacitinib/data/originals/scRNA"
mkdir -p $DATA_DIR

SRA_DATA_DIR="$DATA_DIR/dataset_2/sra"
mkdir -p $SRA_DATA_DIR

OUTPUT_DIR="$DATA_DIR/dataset_2/d2_scRNA_raw_data"
mkdir -p $OUTPUT_DIR

prefetch $SRR --output-directory $SRA_DATA_DIR

fasterq-dump $SRA_DATA_DIR/$SRR/$SRR.sra --split-files --threads 4 -O $OUTPUT_DIR

echo "✅ $SRR download completed!"
