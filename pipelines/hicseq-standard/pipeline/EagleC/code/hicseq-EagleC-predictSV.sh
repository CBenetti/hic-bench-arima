#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH --array=1-16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
##SBATCH --ntasks=1
#SBATCH --job-name=predictSV
#SBATCH --output=predictSV.%A.%a.out
#SBATCH --error=predictSV.%A.%a.err


##
## USAGE: hicseq-template.sh GENOME PREFIX RES1 OUTPUT_FORMAT
##
sleep $((SLURM_ARRAY_TASK_ID * 5))
if [ "$#" -ne 2 ]; then
  grep '^##' "$0"
  exit 1
fi


genome=$1
prefix=$2

module unload python
module load condaenvs/new/EagleC

if cooler info "filtered_cnv_5kb.cool" &>/dev/null && \
   cooler info "filtered_cnv_10kb.cool" &>/dev/null && \
   cooler info "filtered_cnv_50kb.cool" &>/dev/null; then

	predictSV --hic-5k filtered_cnv_5kb.cool \
          --hic-10k filtered_cnv_10kb.cool \
          --hic-50k filtered_cnv_50kb.cool \
          -O $prefix -g $genome --balance-type CNV --output-format full \
          --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

else
	exit 1
fi

module unload condaenvs/new/EagleC
