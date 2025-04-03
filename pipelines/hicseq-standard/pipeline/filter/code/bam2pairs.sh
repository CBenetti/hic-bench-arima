#!/bin/bash

# Load shell settings

##
## USAGE: bam2pairs.sh INPUT-DIR OBJECT OUTDIR GENOME
##

if [ "$#" -ne 4 ]; then
  grep '^##' "$0"
  exit 1
fi

branch=$1
objects=$2
outdir=$3
genome=$4

module unload python
module load anaconda3/gpu/new
conda activate pairtools
bam2pairs -c "inputs/genomes/$genome/chrom_$genome.sizes" "$branch/$objects/feather_output/"$objects"_current/$objects.paired.rmdup.bam" "$outdir/$objects"
conda deactivate
module unload anaconda3
