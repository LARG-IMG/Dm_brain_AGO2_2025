#!/bin/bash
#SBATCH --job-name=bowtie2_map_dm3
#SBATCH --output=sbatch_bwt2_dm3_%j.out
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=18:00:00
#SBATCH --account=<your-account>
#SBATCH --partition=shared

set -uex

mkdir -p map_dm3
mkdir -p bigwigs_dm3
# This script uses trimmed fastq files as the input
# you'll need bowtie2, samtools and deepTools available, do it whichever way you prefer: through module load or conda
# you'll also need to generate bowtie2 index

for fq1 in `ls fq_prepared/*_1.fq.gz`
do
	basen=`basename ${fq1}`
        outn=${basen%_1.fq.gz}
	( bowtie2 -x ~/artem/db/dmel/idx/bwt2/dm5.57 -p 24 -1 ${fq1} -2 ${fq1%_1.fq.gz}_2.fq.gz 2>map_dm3/${outn}_dm3.map.log ) | samtools view -@ 24 -F 4 -b | samtools sort -@ 24 > map_dm3/${outn}_dm3.map.bam
	samtools index map_dm3/${outn}_dm3.map.bam
	bamCoverage --bam map_dm3/${outn}_dm3.map.bam -o bigwigs_dm3/${outn}_all --binSize 50 --normalizeUsing CPM -of bigwig --effectiveGenomeSize 120000000 -p 24
	bamCoverage --bam map_dm3/${outn}_dm3.map.bam -o bigwigs_dm3/${outn}_uq --binSize 50 --minMappingQuality 40 --normalizeUsing CPM -of bigwig --effectiveGenomeSize 120000000 -p 24
	echo "${outn} is done!"
done

echo "All is done!"
