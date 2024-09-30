#!/bin/bash
#SBATCH --job-name=MAPSCinzia # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cinzia.benetti@edu.unito.it # Where to send mail
#SBATCH -J hicbch
#SBATCH --mem=5G
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --output=log_files/MAPSCinzia.log
#SBATCH -p cpu_medium

./run



