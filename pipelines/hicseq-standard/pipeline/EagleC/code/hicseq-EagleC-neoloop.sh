#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
##SBATCH --ntasks=1
#SBATCH --job-name=neoloop
#SBATCH --output=neoloop.log
#SBATCH --error=neoloop.err


module unload python
module load anaconda3/gpu/2023.09
conda activate neoloops

if [ "$#" -ne 1 ]; then
  grep '^##' "$0"
  exit 1
fi

objects=$1
cool_files=$( ls *.cool )

assemble-complexSVs -O $objects -B $objects.combined_neoloop.txt --balance-type CNV --protocol insitu --nproc 4 -H ${cool_files[@]}

neoloop-caller -O $objects.neo-loops.txt --assembly $objects.assemblies.txt --balance-type CNV --protocol insitu --prob 0.95 --nproc 4 -H ${cool_files[@]}
