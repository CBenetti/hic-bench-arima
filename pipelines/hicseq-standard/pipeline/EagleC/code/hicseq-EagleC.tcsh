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

echo $genome
echo $specie
echo $resolution

scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set enzyme = `cut -f1-2,5-7 inputs/sample-sheet.tsv | grep -w "$objects" | cut -f4 | head -n1`
set dir = `pwd`


### run SV prediction
if ( "$resolution" == "5,10,50" ) then
	set branch_10 = `echo "$branch" | sed 's/5/10/g'`
	set branch_50 = `echo "$branch" | sed 's/5/50/g'`
	# create path
	scripts-create-path $outdir/
	ln -ns `realpath "${branch_10}/${objects}/filtered_cnv_10kb.cool"` "$outdir/filtered_cnv_10kb.cool"
	ln -ns `realpath "${branch}/${objects}/filtered_cnv_5kb.cool"` "$outdir/filtered_cnv_5kb.cool"
	ln -ns `realpath "${branch_50}/${objects}/filtered_cnv_50kb.cool"` "$outdir/filtered_cnv_50kb.cool"
	cd $outdir
	sbatch -W "$dir/code/hicseq-EagleC-predictSV.sh" $genome $objects
else
	ln -ns `realpath "${branch}/${objects}/filtered_cnv_${resolution}kb.cool"` "$outdir/filtered_cnv_${resolution}kb.cool"
	cd $outdir
	sbatch -W "$dir/code/hicseq-EagleC-predictSV_singleres.sh" $genome $objects
endif

###create a final results file to use

	$dir/code/SV_merge_results.sh $objects

###annotate gene fusions
if ( "$run_gene_fusion" == "YES" ) then
	$dir/code/hicseq-EagleC-gene_fusion.sh $objects $specie $ens_release
endif


###predict neoloop formation
if ( "$run_neoloop" == "YES" ) then
	sbatch -W $dir/code/hicseq-EagleC-neoloop.sh $objects
endif

cd $dir

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."

