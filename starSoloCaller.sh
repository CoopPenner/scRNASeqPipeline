#!/bin/sh

#  starSoloCaller.sh
#  
#
#  Created by Cooper Penner on 1/13/25.
#  

# now rerunning


FASTQ_DIR="/Volumes/PC60/PDMidBrainFastq_smajic_2022"
GENOME_DIR="/Users/pennerc/Desktop/star_genome_index"
OUTPUT_DIR="/Users/pennerc/Desktop/STARsolo_Output"
WHITELIST="/Users/pennerc/Desktop/refdata-gex-GRCh38-2024-A/whitelist.txt"
GTF_FILE="/Users/pennerc/Desktop/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

mkdir -p "$OUTPUT_DIR"

for sample in $(ls "$FASTQ_DIR" | grep '_R1_001.fastq' | sed 's/_R1_001.fastq//'); do
    echo "Processing sample: $sample"

    # Run STARsolo with correct file order: Read2 (cDNA) first, Read1 (barcode) second
    STAR --runThreadN 8 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$FASTQ_DIR/${sample}_R2_001.fastq" "$FASTQ_DIR/${sample}_R1_001.fastq" \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist "$WHITELIST" \
         --sjdbGTFfile "$GTF_FILE" \
         --soloUMIlen 12 \
         --outFileNamePrefix "$OUTPUT_DIR/${sample}_"
done
