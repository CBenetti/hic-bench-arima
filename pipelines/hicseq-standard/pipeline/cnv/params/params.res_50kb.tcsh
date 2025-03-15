#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = hint  ## hint or neoloop
set resolution = 50    # 50kb resolution. It must match the resolution of one of the cool objects

##If neoloop is chosen, ensure that the cool files have bin generated in the tracks step
##Also, ensure that the files include chr nomenclature
##If not (an error pops up), resolve like this:
## In a bash node:
## module load anaconda3/gpu/2023.09
## conda activate neoloops
## python
###### import pandas as pd
###### import cooler
###### dic={chr:'chr'+str(chr) for chr in c.chromnames}
###### c=cooler.Cooler('matrix.mcool::resolutions/5000') ## replace with file names of interest
###### cooler.rename_chroms(c,dic)
