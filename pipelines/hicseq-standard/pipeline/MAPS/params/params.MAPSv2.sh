#!/bin/bash

# Run the tcsh script and capture the output
Organism="hg38"

# Custom set parameters
bwa__index=inputs/data/genomes/$Organism/genome/bwa.index/genome.fasta
call__peaks=0 # "1" to call peaks using MACS2, "0" to use ChIP peak file provided
peak__type="broad" # Required if call_peaks=1. Type of peaks to call. Must be either "broad" (H3K4me3) or "narrow" (CTCF)
Feather=1 # "1" to run feather, "0" to skip
Maps=1 # "1" to run MAPS on data previously processed with feather, "0" to skip
# setting both feather and MAPS to zero will run the Arima QC metrics only on previous processed data
Threads=16
patterned__flowcell="1" # Use "1" for deep sequencing and "0" for shallow sequencing datasets. "1" if the data was sequenced on a patterned flowcell.
macs2__filepath="inputs/peaks/results/peaks/Arima-MAPS-test/ENCFF247YHM.UW.bed" #Set only if MACS is already run --
                   #-> If so, set to: /inpdir/peaks/results/peaks/${objects}/_peaks.narrowPeak


