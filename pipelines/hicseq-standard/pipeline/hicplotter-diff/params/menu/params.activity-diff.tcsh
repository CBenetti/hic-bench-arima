#!/bin/tcsh

source ./inputs/params/params.tcsh

# this is necessary in order to take care of module conflicts in our system
module unload gcc
module unload python
module load python/2.7.3

# HiCplotter path
set hicplotter_path = ./code/HiCPlotter4.py
set gene_path = `readlink -f $genome_dir/hicplotter-genes.bed`
set chrom_excluded = 'chr[MY]'

#ll $hicplotter_path
#ll $gene_path

# Get the boundary scores and the domains
set branch_short = `echo $branch | sed 's/.*results\///' | sed 's/^matrix-distnorm.[^/]\+\///'`
set group_var = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f2`
set bscores_branch = ../boundary-scores/results/boundary-scores.$group_var.activity_500kb/$branch_short
set domains_branch = ../domains/results/domains.$group_var.hicratio.d_0500/$branch_short
set bdiff_branch = ../boundary-diff/results/boundary-diff.$group_var.standard/$branch_short

set cell_type1 = `./code/code.main/read-sample-sheet.tcsh $sheet "$object1" cell-type`
set cell_type2 = `./code/code.main/read-sample-sheet.tcsh $sheet "$object2" cell-type`
set methods = 'ratio'
set kappas = `cd $bscores_branch/$object1; ls -1 all_scores.k=*.tsv | cut -d'.' -f2`
set bedgraphs = ()
set bedgraph_labels = ()
set beds = ()
set bed_labels = ()

foreach kappa ($kappas)
  set f1 = $bscores_branch/$object1/all_scores.$kappa.tsv
  set f2 = $bscores_branch/$object2/all_scores.$kappa.tsv
  set bfiles = ()
  set blabels = ()
  set dfiles = ()
  set dlabels = ()

  # add bscores for each method
  foreach m ($methods)
    set k = `head -1 $f1 | tr '\t' '\n' | grep -n "^$m"'$' | cut -d':' -f1`             # get column number for method
    cat $f1 | sed '1d' | cut -f1,$k | sed 's/:/\t/' | sed 's/-/\t/' | grep -v '	NA$' >! $workdir/bscores.$object1.$kappa.$m.bedGraph  # convert to bedgraph, remove NA values 
    cat $f2 | sed '1d' | cut -f1,$k | sed 's/:/\t/' | sed 's/-/\t/' | grep -v ' NA$' >! $workdir/bscores.$object2.$kappa.$m.bedGraph  # convert to bedgraph, remove NA values
    set bfiles = ($bfiles bscores.$object1.$kappa.$m.bedGraph bscores.$object2.$kappa.$m.bedGraph)
    set blabels = ($blabels $object1 $object2)
  end

  # add CTCF ChIP-seq for obj #1
  #if (-e inputs/data.external/$cell_type1/CTCF.bedGraph) then
  #  set bfiles = ($bfiles `readlink -f inputs/data.external/$cell_type1/CTCF.bedGraph`)
  #  set blabels = ($blabels $cell_type1-CTCF)
  #endif

  # add CTCF ChIP-seq for obj #2
  #if (-e inputs/data.external/$cell_type2/CTCF.bedGraph) then
  #  set bfiles = ($bfiles `readlink -f inputs/data.external/$cell_type2/CTCF.bedGraph`)
  #  set blabels = ($blabels $cell_type2-CTCF)
  #endif

  # add domains for each one of the samples 
  set d1 = $domains_branch/$object1/domains.$kappa.bed
  set d2 = $domains_branch/$object2/domains.$kappa.bed
  cat $d1 >! $workdir/domains.$object1.$kappa.bed
  cat $d2 >! $workdir/domains.$object2.$kappa.bed
  set dfiles = ($dfiles domains.$object2.$kappa.bed domains.$object1.$kappa.bed)
  set dlabels = ($dlabels $object2 $object1)

  # add to the lists
  set bedgraphs = ($bedgraphs `echo $bfiles | tr ' ' ','`)
  set bedgraph_labels = ($bedgraph_labels `echo $blabels | tr ' ' ','`)
  set beds = ($beds `echo $dfiles`)
  #set bed_labels = ($bed_labels `echo $dlabels | tr ' ' ','`)

end

# regions to plot [TODO: take only $object1.$object2 bdiff?]
set bdiff_files = `find $bdiff_branch/ -name bdiff.k=001.tsv | grep 'T_cell-no_treat' | grep -vE 'GSI'`
set regions = `cat params/diff-activity.bed | cut -f-3 | gtools-regions -center | gtools-regions -shiftp -5p -2000000 -3p +2000000 | gtools-regions -bounds -g $genome_dir/genome.bed | sed 's/\t/:/' | sed 's/\t/-/'`

set tiles = "params/regions.bed"
set tiles_labels = "regions"
set domainbars = 1         # Plot the domains as bars
set highlight = 0
set highlight_bed = ""
set fileheader = 0         # Either 1 or 0 (header / no header)
set insulation_score = 0   # Either 1 or 0 (include insulation index or not)
set compare = 1            # Make the comparisons (based on logFC)
set pair = 0               # Make all the pairwise combinations


