#!/bin/bash

# Load shell settings

##
## USAGE: hicseq-template.sh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if [ "$#" -ne 7 ]; then
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

# Create path

# Save variables

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

mkdir -p "$outdir/neoloop/cnv"

module unload python
module load anaconda3/gpu/2023.09
conda activate neoloops
calculate-cnv -H "$cool_file" -g "$genome" -e "$enzyme" --output "$outdir/neoloop/cnv/$object.$resolution.CNV.bedGraph" --logFile "$outdir/neoloop/cnv/$object.calculate-cnv"
conda deactivate
module unload anaconda3/gpu/2023.09


