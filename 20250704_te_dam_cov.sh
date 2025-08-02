#!/bin/bash

#SBATCH --job-name=damidTEmap
#SBATCH --output=sbatch_damidTE_%j.out
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=06:40:00
#SBATCH --account=<your-account>
#SBATCH --partition=shared

set -uex

# Load modules
ml deepTools/3.5.5
ml hisat2/2.2.1
ml samtools/1.20

# Directories
TE_BAMS_DIR="te_bams"
GENOME_BAMS_DIR="genome_bams"
OUT_DIR="TE_cov_Dam"
mkdir -p "$OUT_DIR"

# Loop over TE BAMs, using the same filenames in genome_bams for normalization
for bam in ${TE_BAMS_DIR}/*DamHP1_Dam.bam; do
  sample=$(basename "$bam")             # e.g., sample_DamHP1.bam
  genome_bam="${GENOME_BAMS_DIR}/${sample}"

  # count mapped reads in the genome BAM (exclude unmapped)
  mapped_reads=$(samtools view -c -F 4 "$genome_bam")

  # compute scale factor for CPM: 1e6 / mapped_reads
  scale=$(awk "BEGIN {printf \"%.6f\", 1000000/$mapped_reads}")

  # Generate bedGraph coverage normalized by CPM matching genome mapping
  bamCoverage \
    -b "$bam" \
    -o "${OUT_DIR}/${sample%.bam}.bg" \
    --scaleFactor $scale \
    --binSize 50 \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK} \
    --outFileFormat bedgraph

  echo "Generated coverage for ${sample%.bam} (CPM normalized using $mapped_reads reads)"
done

