#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compartments.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"


# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# make a list of filtered read files for all input objects
set reg_files = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/filtered.reg.gz"}'`

if ($tool == fithic) then
	./code/hicseq-loops-fithic.tcsh $outdir $params "$reg_files" $genome $branch "$objects"

else if ($tool == fithichip) then
	if (`echo $branch | cut -f 5 -d"/"` == "align.by_sample.hicpro") then 
	      #run PeakInferHiChIP.sh to get peakfile, which is needed by fithichip to call loops 
	      set hicpro = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/hicpro/"}'`
              ./code/hicseq-PeakInferHiChIP-fithichip.tcsh $outdir $params "$hicpro" $genome $branch "$objects"
	
	      #set allvalidpair file
	      set allValidPair = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/hicpro/*allValidPairs"}'`
	      ./code/hicseq-loops-fithichip.tcsh $outdir $params "$allValidPair" $genome $branch "$objects"
	else 
		set branch_short = `echo $branch | cut -d'/' -f4-`
		set files = `find ../tracks/* | grep -c "$branch_short"`
		echo $files
		if(`echo $branch | cut -f 5 -d"/"` == "MAPS.by_sample.MAPSv2" && $files > 0) then
			scripts-create-path "$outdir/PeakInferHiChIP/MACS2_ExtSize/"
			foreach obj ($objects)
				find "../MAPS/results/MAPS.by_sample.MAPSv2/$obj" -type f -name '*Peak' -exec cat {} + >> $outdir/PeakInferHiChIP/MACS2_ExtSize/out_macs2_peaks.narrowPeak
			end
			sed -i 's/chr//g' $outdir/PeakInferHiChIP/MACS2_ExtSize/out_macs2_peaks.narrowPeak
			set hic = `echo $objects | tr ' ' '\n' | awk -v d="../tracks/results/tracks.by_sample.juicer/$branch_short" '{print d"/"$0"/filtered.hic"}'`
			echo $hic
              		./code/hicseq-loops-fithichip.tcsh $outdir $params "$hic" $genome $branch "$objects"
			awk '{if (NR>1) {if (substr($1,1,1) ~ /^[0-9]/ ) {print "chr"$1"\t"$2"\t"$3"\tchr"$4"\t"$5"\t"$6"\t"$7} else {print $0}}}' $outdir/FitHiChIP/FitHiChIP_Peak2ALL_b"$winsize"_L"$mindist"_U"$maxdist"/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q"$qval"_MergeNearContacts.bed > $outdir/loops_filtered_bias_raw.bedpe
			awk 'NR>1' ${outdir}/FitHiChIP/FitHiChIP_Peak2ALL_b"$winsize"_L"$mindist"_U"$maxdist"/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q"$qval"_MergeNearContacts.bed | cut -f 9 >! ${outdir}/qval.txt
			cp $outdir/loops_filtered_bias_raw.bedpe $outdir/loops_filtered_nobias_raw.bedpe
			cat $reg_files | gunzip >! $outdir/filtered.reg
			set intra_reads = `cat $outdir/filtered.reg | awk '$2 == $6' | wc -l`
			awk -v var="$intra_reads" '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/(var/1000000) }' $outdir/loops_filtered_nobias_raw.bedpe >! $outdir/loops_filtered_nobias_cpm.bedpe
			cp $outdir/loops_filtered_nobias_cpm.bedpe $outdir/loops_filtered_bias_cpm.bedpe
			paste ${outdir}/loops_filtered_nobias_cpm.bedpe ${outdir}/qval.txt >! ${outdir}/loops_labeled_qval.bedpe
			cp $outdir/FitHiChIP/Summary_results_FitHiChIP.html $outdir/
		else
	      		echo "Error: Fithichip loop calling requires hic-pro or MAPS output." | scripts-send2err
		endif
	endif

else
  	echo "Error: Loops calling tool $tool not supported." | scripts-send2err
endif


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
