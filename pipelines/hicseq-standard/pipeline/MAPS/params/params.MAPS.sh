#!/bin/bash

# Source the params file
source ./inputs/params/params.sh

# Run the read-sample-sheet script and assign the output to a variable
Organism=$(./code/read-sample-sheet.sh "$sheet" "$object" genome)

# Set other variables
bwa__index="inputs/data/genomes/$genome/genome/bwa.index/genome.fasta"
macs2__filepath="inpdirs/peaks/results/peaks.by.sample"
fastq__dir="inputs/fastq"
bin__size=10000
binning__range=1000000
length__cutoff=1000
Threads=16
Models="pospoisson"  # or "negbinom"
sex__chroms__to__process="X"  # or "Y", "XY", "NA"
generate__hic=1  # or 0

# Get the paths to Python and Rscript
python__path=$(which python)
Rscript__path=$(which Rscript)
