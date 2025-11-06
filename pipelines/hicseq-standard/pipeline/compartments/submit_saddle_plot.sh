#!/bin/bash
#SBATCH -J saddle_plot
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH -N 1

##usage: select two branches (one for hic and one for compartments) and modify saddle_plots.R
##	then chose a branch and modify the following line to list the samples

branch = "results/compartments.by_sample.homer.res_100kb/filter.by_sample.mapq_20_mindist0/align.by_sample.hicpro"

for sample in "$branch"/*; do
  sample_name=$(basename "$sample")
  sbatch --export=SAMPLE="$sample_name" saddle_plots.sh
done
