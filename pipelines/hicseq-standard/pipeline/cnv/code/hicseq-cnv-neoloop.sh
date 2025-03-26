#!/bin/bash

# Load shell settings

##
## USAGE: hicseq-template.sh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if [ "$#" -ne 8 ]; then
  grep '^##' "$0"
  exit 1
fi

outdir="$1"
cool_file="$2"
object="$3"
enzyme="$4"
genome="$5"
params="$6"
resolution="$7"
chr_rename="$8"

# Create path

# Save variables

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

mkdir -p "$outdir/neoloop/cnv"

module unload python
module load anaconda3/gpu/2023.09
conda activate neoloops
if [[ "$cool_file" == *.mcool ]]; then
	resolution_bp=$((resolution * 1000))
	cool_file="$cool_file::resolutions/$resolution_bp"
fi

if cooler info "$cool_file" &>/dev/null; then
	if [[ "$chr_rename" == "YES"]]; then
		python ./code/convert_dict.py $cool_file
	fi 
	calculate-cnv -H "$cool_file" -g "$genome" -e "$enzyme" --output "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --logFile "$outdir/neoloop/cnv/$object.calculate-cnv"
	cp $cool_file $outdir/filtered_cnv_"$resolution"kb.cool
	correct-cnv -H "$outdir/filtered_cnv_"$resolution"kb.cool" --cnv-file "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --nproc 4 -f
else
	echo "Resolution group does not exist: $cool_file"
fi	
conda deactivate
module unload anaconda3/gpu/2023.09


