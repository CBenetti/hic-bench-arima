#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-cnv.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

set res = `echo $params | sed -E 's/.*res_([0-9]+kb).*/\1/'`

echo $outdir
echo $params
echo $objects
echo $branch

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set hic_file = $branch/$objects/filtered.hic
set enzyme = `cut -f1-2,5-7 inputs/sample-sheet.tsv | grep -w "$objects" | cut -f4 | head -n1`
set cool_file = $branch/$objects/filtered.cool
if ($tool == hint) then
  ./code/hicseq-cnv-hint.tcsh $outdir $params $genome $enzyme $hic_file
else
   if ($tool == neoloop)then
      ./code/hicseq-cnv-neoloop.sh $outdir $cool_file $objects $enzyme $genome $params $resolution
   else
     echo "Error: Cnv tool $tool not supported." | scripts-send2err
   endif
endif


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
