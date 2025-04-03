#--------------------------------
#  USAGE: 1) modify gene name line (line 43)
#         2) module load condaenvs/new/neoloop
#         3) python3 code/plot.py 'OUTDIR' 'SAMPLE' 'RESOLUTION(Kb)' ASSEMBLY-LINE-TO-PLOT 'FIGURE NAME'
#--------------------------------

##########################################
### 1 - IMPORT PACKAGES AND READ INPUT VAR
##########################################
import sys
import os
import pandas as pd
import cooler
from neoloop.visualize.core import *
from pathlib import Path
sample = sys.argv[2]
res = sys.argv[3]
figname = sys.argv[5]
directory = sys.argv[1]
assembly_line = sys.argv[4]
cool_f = os.path.join(directory, 'filtered_cnv_'+res+'kb.cool')
clr = cooler.Cooler(cool_f)
assembly = os.path.join(directory,sample+'.assemblies.txt')
with open(assembly, "r") as f:
    lines = f.readlines()  # Read all lines into a list
    l = lines[int(assembly_line)].strip()

outfig = os.path.join(directory,figname)
directory = Path(directory)
cnv = os.path.join('inpdirs','cnv','results',Path(*directory.parts[2:]),'neoloop','cnv',sample+'.'+res+'.CNV.bw')
neoloop = os.path.join(directory,sample+'.neo-loops.txt')

##############
### 2 - PLOT
##############

vis = Triangle(clr, l, n_rows=4, figsize=(7, 5.3), track_partition=[5, 0.4, 0.4, 0.5], correct='weight', span=1400000, space=0.03)
vis.matrix_plot(vmin=0, cbr_fontsize=9)
vis.plot_chromosome_bounds(linewidth=2)
vis.plot_loops(neoloop, face_color='none', marker_size=40, cluster=True, onlyneo=True)

########## LINE TO MODIFY HERE ! ############
vis.plot_genes(filter_=['CAAP1','DENND4C','SLC2A2','ACER2'],label_aligns={'DENND4C':'right'},fontsize=10)
#############################################

vis.plot_signal('CNV', cnv, label_size=10, data_range_size=9, max_value=1.8, color='#6A3D9A')
vis.plot_chromosome_bar(name_size=13, coord_size=10, color_by_order=['#1F78B4','#33A02C'])
vis.outfig(outfig, dpi=800)
