#!/bin/bash
#SBATCH -J resolution
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH -N 1

module load anaconda3/gpu/2023.09
conda activate neoloops

##usage: chose resolutions to calculate and modify calculate_resolution.py
##	then chose a branch and modify the following line

branch = "results/tracks.by_sample.juicer_cool.multires/filter.by_sample.mapq_20_mindist0/align.by_sample.hicpro"

for sample in "$branch"/*; do
  sample_name=$(basename "$sample")
  python3 calculate_resolution.py "$sample_name" "$branch"
done
