#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH --array=1-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
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
file=$( ls *.cool )

module unload python
module load condaenvs/new/EagleC
if cooler info "${file[0]}" &>/dev/null; then
	predictSV-single-resolution --hic $file -O "$prefix.combined.txt" -g $genome --balance-type CNV --output-format full
else
	echo "no such file or directory"
	exit 1
fi
module unload condaenvs/new/EagleC
