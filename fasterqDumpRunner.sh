#!/bin/sh

#  fasterqDumpRunner.sh
#  
#
#  Created by Cooper Penner on 1/13/25.
#  
# Define input file and output directory
OUTPUT_DIR="/Volumes/PC60/PDMidBrainFastq_smajic_2022"
TEMP_DIR="${OUTPUT_DIR}/temp"
SRR_LIST="${OUTPUT_DIR}/SRR_Acc_List.txt"

# Check if the SRR list file exists
if [ ! -f "$SRR_LIST" ]; then
    echo "Error: SRR list file not found at $SRR_LIST"
    exit 1
fi

# Create the output and temporary directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# Loop through each SRR ID and run fasterq-dump
while IFS= read -r SRR_ID; do
    echo "Processing $SRR_ID..."
    fasterq-dump "$SRR_ID" --threads 8 --outdir "$OUTPUT_DIR" --temp "$TEMP_DIR"
    echo "$SRR_ID completed."
done < "$SRR_LIST"

echo "All SRR IDs have been processed!"
