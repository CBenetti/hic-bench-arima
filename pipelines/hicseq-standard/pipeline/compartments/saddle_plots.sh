#!/bin/bash
#SBATCH -J saddle_plots
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH -N 1


module load r/4.3.2
Rscript saddle_plots.R "$SAMPLE"
