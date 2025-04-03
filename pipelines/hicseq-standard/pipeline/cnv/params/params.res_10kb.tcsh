#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = neoloop  ## hint or neoloop
set resolution = 10    # 10kb resolution. It must match the resolution of one of the cool objects, and be expressed in kb
set chr_rename = "YES"  #Set to YES ONLY if neoloop method is selected and ONLY if it's the first time running the module on the dataset. DO NOT LEAVE EMPTY
