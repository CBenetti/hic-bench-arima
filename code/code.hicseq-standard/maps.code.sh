#!/bin/bash

# Assign command-line arguments to variables
outd="$1"
params="$2"
branch="$3"
# Convert space-separated arguments to an array
IFS=' ' read -r -a objects <<< "$4"
# Set all possible customizable parameters

source "$params"
echo "$Organism"
echo "$outd"
echo "./inputs/fastq/${objects[@]}/${objects[@]}"
fastq1=$(./code/read-sample-sheet.tcsh inputs/sample-sheet.tsv "${objects[@]}" fastq-r1 | tr ',' '\n' | sed "s|[^ ]*|"inputs/fastq/"&|g")
zcat $fastq1 | gzip > inputs/fastq/"${objects[@]}"/"${objects[@]}"_R1.fastq.gz
fastq2=$(./code/read-sample-sheet.tcsh inputs/sample-sheet.tsv "${objects[@]}" fastq-r2 | tr ',' '\n' | sed "s|[^ ]*|"inputs/fastq/"&|g")
zcat $fastq2 | gzip > inputs/fastq/"${objects[@]}"/"${objects[@]}"_R2.fastq.gz

./bin/Arima-MAPS_v2.0.sh -C ${call__peaks} -p ${peak__type} -F ${Feather} -M ${Maps} \
-I "./inputs/fastq/${objects[@]}/${objects[@]}" -O ${outd} -m ${macs2__filepath} \
-o ${Organism} -b ${bwa__index} -t ${Threads} -f ${patterned__flowcell}
