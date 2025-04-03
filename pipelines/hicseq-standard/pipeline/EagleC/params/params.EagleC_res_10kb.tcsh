#!/bin/tcsh

source ./inputs/params/params.tcsh

set genome = "hg38" ## hg38,hg19,mm...
set resolution = 10 # 10kb resolution. It must match the resolution of one of the cool objects, and be expressed in kb.IF using a single resolution different from 10, also modify run script ! For muliple resolution, provide ONLY the following comma separated list: 5,10,50
set run_gene_fusion = "YES" #Set to YES ONLY if gene fusions needed
set run_neoloop = "YES" #Set to YES ONLY if neoloop  needed
set specie = "human"
set ens_release = 110 #110 if hg38, 73 if hg19
