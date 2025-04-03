#!/bin/bash
module unload python
module load condaenvs/new/EagleC

objects=$1
txt_files=( *.txt )
if [ -e "${txt_files[0]}" ]; then
	final=( *combined* )
	if [ -e "${final[0]}" ]; then
		cp "${final[0]}" $objects.combined.txt
		merge-redundant-SVs  --full-sv-files ${final[@]} --output-format NeoLoopFinder -O $objects.combined_neoloop.txt
	else
		target=$objects.combined_neoloop.txt
		merge-redundant-SVs --full-sv-files ${txt_files[@]} --output-format NeoLoopFinder -O $target
		target=$objects.combined.txt
		merge-redundant-SVs --full-sv-files ${txt_files[@]} --output-format full -O $target
	fi
else
	echo "An error occurred"
	exit 1
fi
