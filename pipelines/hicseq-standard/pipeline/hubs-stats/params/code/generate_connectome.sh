#!/bin/bash
module load r/4.4.1
module load bedops/2.4.41
Rscript code/connectome.R 1000
sort-bed peaks.bed > peaks_tmp.bed
awk -v OFS="\t" '{print $1,$2,$3,0,$4}' peaks_tmp.bed > peaks_sorted.bed
awk -v OFS="\t" '{print $1,$3-500,$3-499,$5}' peaks_sorted.bed > viewpoints.bed
module unload r/4.4.1
module unload bedops/2.4.41
awk 'BEGIN{OFS="\t"}{ split($4, a, "_"); $4 = a[2]"_"a[1]"_"a[3]; if ($6 == "+") { start = $2-5000; end = $2+500;} else if ($6 == "-") { start = $2-500; end   = $2+5000;} if (start < 0) start = 0; print $1, start, end, $4, $5, $6;}' ../inputs/genomes/hg38/protein-coding-tss.bed > tss.bed #adjust according to preference
