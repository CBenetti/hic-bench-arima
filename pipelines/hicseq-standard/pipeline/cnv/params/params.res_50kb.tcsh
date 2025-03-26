#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = hint  ## hint or neoloop
set resolution = 50    # 50kb resolution. It must match the resolution of one of the cool objects, and be expressed in kb
set chr_rename = ""  #Set to YES ONLY if neoloop method is selected and ONLY if it's the first time running the module on the dataset
