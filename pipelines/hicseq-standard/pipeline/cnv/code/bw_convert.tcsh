#!/bin/tcsh
#SBATCH -J CNV_bw_convert
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH -N 1


source ./code/code.main/custom-tcshrc

##
## USAGE: bw_convert.tcsh GENOME
##
module load ucscutils/398
set bedG_files = `find ./results -type f -name "*CNV.bedGraph"`
set genome = ./inputs/genomes/$1/chrom_$1.sizes
foreach file ($bedG_files)
	set dir = `dirname $file`
	sort -k1,1 -k2,2n $file > $dir/sorted.bedGraph
	set out = `echo "$file" | sed 's/CNV.bedGraph/CNV.bw/g'`
	echo $genome
	echo $out
	bedGraphToBigWig $dir/sorted.bedGraph $genome $out
end
