#!/bin/bash
#SBATCH --job-name=rename_dataset1_fastqs
#SBATCH --array=1-4
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/prj/bmb_tofacitinib/research_project/scRNA_alignment/scripts/slurms/rename_fastq/dataset1/rename_fastq-%A_%a.out

BASE_DATA_DIR="/scratch/prj/bmb_tofacitinib/data"

# Data where dataset1 fastq lie
DATA_DIR="${BASE_DATA_DIR}/originals/scRNA/dataset_1/d1_scRNA_raw_data"

# Directory where renamed fastq files will be outputed
OUTPUT_DIR="${BASE_DATA_DIR}/originals/scRNA/dataset_1/renamed_fastq"
mkdir -p $OUTPUT_DIR

# Metadata tsv directory
METADATA="/scratch/prj/bmb_tofacitinib/research_project/scRNA/originals/public_metadata.tsv"

# Extract each line from the metadata file containing the sample name and associated SRR files based on the array job iteration 
# 4 array jobs means 4 GSM samples
# +1 to account for the headers

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")

# Split the tsv data into field 1 GSM ID and field 2 SRR IDs
GSM_ID=$(echo "$LINE" | cut -f1)
SRR_LIST=$(echo "$LINE" | cut -f2)

echo "Printing field 1 and 2... (GSM and SRRs)"
echo "$GSM_ID"
echo "$SRR_LIST"

IFS=',' read -ra SRRS <<< "$SRR_LIST"

lane_num=1

for SRR in "${SRRS[@]}"; do
    echo "  Handling SRR: $SRR (Lane $lane_num)"

    R1_SRC="${DATA_DIR}/${SRR}_1.fastq"
    R2_SRC="${DATA_DIR}/${SRR}_2.fastq"

    R1_DEST="${OUTPUT_DIR}/${GSM_ID}_S1_L$(printf "%03d" $lane_num)_R1_001.fastq.gz"
    R2_DEST="${OUTPUT_DIR}/${GSM_ID}_S1_L$(printf "%03d" $lane_num)_R2_001.fastq.gz"

    if [[ -f "$R1_SRC" && -f "$R2_SRC" ]]; then
        gzip -c "$R1_SRC" > "$R1_DEST"
        gzip -c "$R2_SRC" > "$R2_DEST"
        echo "    Renamed and compressed lane $lane_num."
    else
        echo "    WARNING: Missing $R1_SRC or $R2_SRC!"
    fi

    lane_num=$((lane_num + 1))
done

echo "Sample $GSM_ID renaming completed."
