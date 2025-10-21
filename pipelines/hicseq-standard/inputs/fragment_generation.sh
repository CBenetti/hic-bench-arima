#!/bin/bash
#SBATCH --job-name=digest_Arima
#SBATCH --output=digest_Arima_%j.log
#SBATCH --error=digest_Arima_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=01:00:00

  module unload r
  module load r/4.0.3
  module unload python
  module load python/cpu/3.7.2
  module load hic-pro

# Paths
HICPRO_PATH=/gpfs/share/apps/hic-pro/3.0.0/HiC-Pro_3.0.0/bin/utils
GENOME_FASTA=genomes/hg38/bowtie2.index/genome.fa
OUT_BED=Arima_4e_fragments.bed

# Run digestion
$HICPRO_PATH/digest_genome.py -r ^GATC,G^ANTC,C^TNAG,T^TAA  -o $OUT_BED $GENOME_FASTA
