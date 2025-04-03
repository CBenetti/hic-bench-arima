#!/bin/tcsh

module load r/4.2.2
if (! -e ../align/results/.db) then
                mkdir ../align/results
                cp -r results/.db/ ../align/results/
        else
            	tail -n +2 results/.db/db.tsv >> ../align/results/.db/db.tsv
                cat results/.db/run >> ../align/results/.db/run
                #cat results/.db/run.outofdate >> ../align/results/.db/run.outofdate
                Rscript code/combine_rdata.R results/.db/db.RData ../align/results/.db/db.RData ../align/results/.db/db.RData
        endif
	ln -ns ../../MAPS/results/MAPS.by_sample.MAPSv2 ../align/results/MAPS.by_sample.MAPSv2

