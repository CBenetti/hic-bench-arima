#!/bin/tcsh
#SBATCH --job-name=MAPSCinzia # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cinzia.benetti@edu.unito.it # Where to send mail
#SBATCH -J hicbch
#SBATCH --mem=5G
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --output=log_files/MAPSCinzia.log
#SBATCH -p cpu_medium

#Before running this: change parameters and run the pipeline as normal, to ensure that the run.outofdate file is created and executable
#To save previous arcplots and metaplots, MOVE THEM!
code/code.main/scripts-submit-jobs "./results/.db/run.outofdate" 3
