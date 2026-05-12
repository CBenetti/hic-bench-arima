#!/bin/tcsh

# unload tools that may cause conflicts with other tools
module unload r
module unload gcc
module unload samtools
module unload java
module unload python

# load all tools
module load gtools/3.0.0
#module load kentutils/329
module load ucscutils/368
module load samtools/1.20-new
module load bedtools/2.30.0
module load picard-tools/1.88
module load bowtie2/2.5.3
module load r/4.4.2

samtools --version >& /dev/null
if ( $status != 0 ) then
    module unload samtools/1.20-new
    module load samtools/1.9-new
endif

#module load java/9.0.4
# module load macs/2.0.10.20131216

# set sample sheet
set sheet = inputs/sample-sheet.tsv

