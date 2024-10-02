#!/bin/bash
#SBATCH -A b1042
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=predictSV
#SBATCH --error=predictSV.%j.%N.err

module load condaenvs/new/EagleC
predictSV --hic-5k $outdir/matrix_5000.cool \
          --hic-10k $outdir/matrix_10000.cool \
          --hic-50k $outdir/matrix_50000.cool \
          -O $prefix -g $genome --balance-type CNV --output-format full \
          --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

module unload condaenvs/new/EagleC
