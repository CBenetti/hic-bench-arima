#!/bin/bash
module unload python
module load anaconda3/gpu/2023.09
conda activate EagleC

if [ "$#" -ne 3 ]; then
  grep '^##' "$0"
  exit 1
fi

objects=$1
specie=$2
ens_release=$3

if [ -e "$objects.combined.txt" ]; then
	# Substitutes 'chr' with empty string in all lines except the first and saves it in a new file
	awk 'NR==1 {print; next} {gsub(/chr/, ""); print}' "$objects.combined.txt" > "$objects.combined_chr.txt"
	annotate-gene-fusion --sv-file $objects.combined_chr.txt \
                       --output-file $objects.gene-fusions.txt \
                       --buff-size 10000 --skip-rows 1 --ensembl-release $ens_release --species $specie
else

	echo "something went wrong before annotate gene fusion step"
	exit 1
fi
