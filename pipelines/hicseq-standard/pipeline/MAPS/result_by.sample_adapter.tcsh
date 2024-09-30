#!/bin/tcsh
#SBATCH --job-name=MAPSadapterCinzia # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cinzia.benetti@edu.unito.it # Where to send mail
#SBATCH -J MAPS
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --output=log_files/MAPSadapterCinzia.log

source ./code/code.main/custom-tcshrc      # customize shell environment
set winsize = 5 ##In kb, as suggested by MAPS
set mindist = 1000 ##as suggested by MAPS, general parameter
set maxdist = 2000000 ##as suggested by MAPS for loop calling

#adapt output for subsequent analysis

set out_path = "feather/MAPS.by_sample.MAPSv2"
foreach obj (results/MAPS.by_sample.MAPSv2/*/)
  set objects = `basename $obj`
  scripts-create-path $out_path/$objects/
  set current_dir = `pwd`
  set parent_dir = `dirname $current_dir`
  awk '{ \
    $2 = ($2 == "0") ? "+" : ($2 == "16") ? "-" : $2 \
    $6 = ($6 == "0") ? "+" : ($6 == "16") ? "-" : $6 \
    gsub(":", ".", $1); print $1"\t"$3" "$2" "$4" "$4" "$7" "$6" "$8" " $8 \
  }' results/MAPS.by_sample.MAPSv2/$objects/feather_output/"$objects"_current/$objects.hic.input  | sort -k2,2 -k4,4n -k6,6 -k8,8n > $out_path/$objects/filtered.reg
  gzip --force $out_path/$objects/filtered.reg

  if (-e results/MAPS.by_sample.MAPSv2/$objects/MAPS_output) then
    set out_path = "maps/loops.by_sample.MAPS.res_5kb/MAPS.by_sample.MAPSv2"
    scripts-create-path $out_path/$objects/
##Adapt files to the following format
    cp results/MAPS.by_sample.MAPSv2/$objects/MAPS_output/"$objects"_current/"$objects"."$winsize"k.2.sig3Dinteractions.bedpe $out_path/$objects/
    awk -v winsize=$winsize '{ \
      if ( NR == 1 ) \
                print "chr1\tfragmentMid1\tchr2\tfragmentMid2\tcontactCount\tp-value\tq-value\tbias1\tbias2" \
      else \
          	print $1"\t"$2+winsize*500-1"\t"$4"\t"$5+winsize*500-1"\t"$7"\t"sprintf("%e", $9)"\t"sprintf("%e", $9)"\t"sprintf("%e", 1)"\t"sprintf("%e", 1) \
    }' $out_path/$objects/"$objects"."$winsize"k.2.sig3Dinteractions.bedpe >! $out_path/$objects/loops_unfiltered_nobias_raw.tsv

###Calculate normalization
    set intra_reads = `less feather/MAPS.by_sample.MAPSv2/$objects/filtered.reg.gz | awk '$2 == $6' | wc -l`
###Create CPM normalized unfiltered loops files (no bias)
    awk -v var="$intra_reads" '{                                                              \
        if ( NR == 1 )                                                                    \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9                  \
        else                                                                              \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9    \
    }' $out_path/$objects/loops_unfiltered_nobias_raw.tsv | gzip >! $out_path/$objects/loops_unfiltered_nobias_cpm.tsv.gz

##Create filtered and normalized file
  awk -F "\t" -v min="$mindist" -v max="$maxdist" 'NR == 1 || \!/fragment/ && ($4+1-$2) < max && ($4+1-$2) > min' $out_path/$objects/loops_unfiltered_nobias_raw.tsv >! $out_path/$objects/loops_filtered_nobias_raw.tsv
  awk -v intra="$intra_reads" -v var="$winsize" '{ if (NR > 1) print $1"\t"($2+1-var*500)"\t"($2+var*500)"\t"$3"\t"($4+1-var*500)"\t"($4+var*500)"\t"$5/(intra/1000000)}' $out_path/$objects/loops_filtered_nobias_raw.tsv >! $out_path/$objects/loops_filtered_nobias_cpm.bedpe
  endif
end
cd $parent_dir
ln -ns ../MAPS/feather align/results
cp -r MAPS/results/.db/ MAPS/feather/
ln -ns ../MAPS/maps loops/results
cp -r MAPS/results/.db/ MAPS/maps/

