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

resolution_bp=$((resolution * 1000))

if [[ "$cool_file" == *.mcool ]]; then
	cool_file="$cool_file::resolutions/$resolution_bp"
	cooler cp $cool_file $outdir/filtered_cnv_"$resolution"kb.cool
else
	cp $cool_file $outdir/filtered_cnv_"$resolution"kb.cool
fi

if (cooler info "$outdir/filtered_cnv_"$resolution"kb.cool" &>/dev/null); then
	if [[ "$chr_rename" == "YES" ]]; then
		python ./code/convert_dict.py "$outdir/filtered_cnv_"$resolution"kb.cool"
	fi
	calculate-cnv -H "$outdir/filtered_cnv_"$resolution"kb.cool" -g "$genome" -e "$enzyme" --output "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --logFile "$outdir/neoloop/cnv/$object.calculate-cnv"
	segment-cnv --cnv-file "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --binsize $resolution_bp --ploidy 2 --output "$outdir/neoloop/cnv/$object.$resolution.CNV-seg.bedGraph" \
	--nproc 4 --logFile "$outdir/segment_cnv.log"
	plot-cnv --cnv-profile "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --cnv-segment "$outdir/neoloop/cnv/$object.$resolution.CNV-seg.bedGraph" \
	--output-figure-name "$outdir/$object.$resolution.kb.CNV.genome-wide.png" --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 --label-size 7 --tick-label-size 6 --clean-mode
	cooler balance "$outdir/filtered_cnv_"$resolution"kb.cool"
	correct-cnv -H "$outdir/filtered_cnv_"$resolution"kb.cool" --cnv-file "$outdir/neoloop/cnv/$object.$resolution.CNV-seg.bedGraph" --nproc 4 --logFile "$outdir/correct_cnv.log" -f
else
	echo "Resolution group does not exist: $cool_file"
fi
conda deactivate
module unload anaconda3/gpu/2023.09


