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

# === CONFIGURATION ===
SAMPLES="damhp1_samples.txt"  
# pregenerated hisat2 dm6 index location
GENOME_IDX="~/db/dmel/idx/hs2/bdgp6.32"
# pregenerated hisat2 index on TE consensus sequences location
TE_IDX="~/db/dmel/idx/hs2/TE_cons"
THREADS=12

mkdir -p genome_bams te_bams ratios

while read -r damhp1_r1 dam_r1; do
  damhp1_r2=${damhp1_r1/_1.fq.gz/_2.fq.gz}
  dam_r2=${dam_r1/_1.fq.gz/_2.fq.gz}
  prefix=$(basename "$damhp1_r1" _1.fq.gz)

  echo ">>> $prefix"

  # 1) map to genome
  gbam_hp1=genome_bams/${prefix}_DamHP1.bam
  if [[ -f $gbam_hp1 ]]; then
	  echo " [skip] ${gbam_hp1} already exists!"
  else
  hisat2 -p $THREADS -x $GENOME_IDX \
      -1 "$damhp1_r1" -2 "$damhp1_r2" \
    | samtools view -bS -@2 > genome_bams/${prefix}_DamHP1.bam
  fi
  
  gbam_dam=genome_bams/${prefix}_Dam.bam
  if [[ -f $gbam_dam ]]; then
	  echo " [skip] ${gbam_dam} already exists!"
  else
    hisat2 -p $THREADS -x $GENOME_IDX \
        -1 "$dam_r1" -2 "$dam_r2" \
      | samtools view -bS -@2 > genome_bams/${prefix}_Dam.bam
  fi

  # 2) count mapped reads via samtools view
  mapped_hp1=$(samtools view -c -F 4 genome_bams/${prefix}_DamHP1.bam)
  mapped_dam=$(samtools view -c -F 4 genome_bams/${prefix}_Dam.bam)
  echo "  Mapped reads Dam-HP1: $mapped_hp1, Dam-only: $mapped_dam"

  # 3) map to TE consensus
  for sample in DamHP1 Dam; do
    fq1_var="${sample,,}_r1"
    fq2_var="${sample,,}_r2"
    tebam=te_bams/${prefix}_${sample}.bam
    if [[ -f ${tebam} ]]; then
	    echo " [skip] ${tebam} already exists!"
    else
      hisat2 -p $THREADS -x $TE_IDX \
        -1 "${!fq1_var}" -2 "${!fq2_var}" \
      | samtools view -bS -@8 -F 4 - \
      | samtools sort -@8 - > te_bams/${prefix}_${sample}.bam
      echo "  TE mapping done for $sample"
    fi
  done

  # 3.1) index TE bams
  for sample in DamHP1 Dam; do
	  tebam=te_bams/${prefix}_${sample}.bam
	  if [[ -f ${tebam}.bai ]]; then
		  echo " [skip] ${tebam} is already indexed!"
	  else
		  samtools index ${tebam}
	  fi
  done

  # 4) compute scale factor and bamCompare with log2 ratio
  scale_hp1=$(awk "BEGIN {printf \"%.6f\", $mapped_dam / $mapped_hp1}")
  bgout=ratios/${prefix}_log2Dam.bedgraph
  if [[ -f ${bgout} ]]; then
	  echo " [skip] ${bgout} already exists!"
  else
    bamCompare \
      -b1 te_bams/${prefix}_DamHP1.bam \
      -b2 te_bams/${prefix}_Dam.bam \
      --scaleFactorsMethod None \
      --scaleFactors $scale_hp1:1 \
      --operation log2 \
      --outFileFormat bedgraph \
      -o ratios/${prefix}_TE_log2Dam.bedgraph
  fi

  echo "  log2ratio track: ratios/${prefix}_TE_log2Dam.bedgraph"
done < "$SAMPLES"
