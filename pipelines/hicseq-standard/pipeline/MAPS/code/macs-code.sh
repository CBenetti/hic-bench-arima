#!/bin/bash

# Assign command-line arguments to variables
outd="$1"
params="$2"
branch="$3"
# Convert space-separated arguments to an array
IFS=' ' read -r -a objects <<< "$4"
# Set all possible customizable parameters

source "$params"

Arima-MAPS_v2.0.sh -C ${call__peaks} -p ${peak__type -F ${Feather} -M ${Maps}\
 -I ./inputs/fastq/${objects[@]} -O "${outd}/$objects[@]"\
 -m ${macs2__filepath} -o ${Organism} -b ${bwa__index} -t ${Threads} -f ${patterned__flowcell}
