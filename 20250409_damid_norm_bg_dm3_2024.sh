#!/bin/bash

#SBATCH --job-name=bamcomparetest
#SBATCH --output=sbatch_bct_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:30:00
#SBATCH --account=<your-account>
#SBATCH --partition=shared

# 2024-09-09 
# A script to normalize DamID bams getting bigwig files in the output
set -uex
# Define the directory where BAM files are located
BAM_DIR="map_dm3"
OUT_DIR="norm_bg_dm3"  # Directory for output files
GENOME_SIZE=120000000        # Effective genome size
BIN_SIZE=50                  # Bin size
NUM_THREADS=8                # Number of processors

# Create the output directory if it doesn't exist
mkdir -p $OUT_DIR

# Loop over all BAM files that match the pattern *DamHP1*.bam
for bam_file in $BAM_DIR/*DamHP1*.bam; do
    # Extract the base name of the file without extension
    base_name=$(basename "$bam_file" .bam)

    # Extract the prefix and replicate number (assuming the format like E1C_DamHP1_rep1_dm3.map.bam)
    prefix=$(echo $base_name | cut -d'_' -f1)  # This extracts the first part, e.g., E1C
    replicate=$(echo $base_name | cut -d'_' -f3)  # This extracts the replicate number, e.g., rep1

    # Construct the corresponding Dam file name
    dam_file="${BAM_DIR}/${prefix}_Dam_${replicate}_dm3.map.bam"

    # Check if the corresponding Dam BAM file exists
    if [ -f "$dam_file" ]; then
        echo "Normalizing $bam_file to $dam_file..."

        # Define output file name
        output_file_uniq="${OUT_DIR}/${base_name}_vs_${prefix}_Dam_${replicate}_log2_CPM_uniq.bedgraph"
	output_file_multi="${OUT_DIR}/${base_name}_vs_${prefix}_Dam_${replicate}_log2_CPM_multi.bedgraph"

        # Run bamCompare for normalization
        bamCompare \
            --bamfile1 "$bam_file" \
            --bamfile2 "$dam_file" \
            --outFileName "$output_file_uniq" \
            --outFileFormat bedgraph \
            --scaleFactorsMethod None \
            --operation log2 \
            --normalizeUsing CPM \
            --effectiveGenomeSize $GENOME_SIZE \
	    --minMappingQuality 42 \
            --binSize $BIN_SIZE \
            --numberOfProcessors $NUM_THREADS

        echo "Output saved to $output_file_multi"

        bamCompare \
            --bamfile1 "$bam_file" \
            --bamfile2 "$dam_file" \
            --outFileName "$output_file_multi" \
            --outFileFormat bedgraph \
            --scaleFactorsMethod None \
            --operation log2 \
            --normalizeUsing CPM \
            --effectiveGenomeSize $GENOME_SIZE \
            --minMappingQuality 0 \
            --binSize $BIN_SIZE \
            --numberOfProcessors $NUM_THREADS

        echo "Output saved to $output_file_uniq"
    else
        echo "Dam file for $bam_file not found: $dam_file"
    fi
done

echo "All files processed."

