#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=sbatch_trim_galore%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=<your-account>
#SBATCH --partition=shared

set -uex

# Output directory (adjust if necessary)
OUTDIR="trimmed_fq"

# Create the output directory if it doesn't exist
mkdir -p $OUTDIR

for fq in fastq/*_1.fq.gz
do
	trim_galore --paired --cores 8 --output_dir $OUTDIR ${fq} ${fq%_1.fq.gz}_2.fq.gz
done

echo "Trimming completed" 

