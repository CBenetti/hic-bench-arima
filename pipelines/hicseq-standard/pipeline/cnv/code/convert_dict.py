import sys
import pandas as pd
import cooler
filename = sys.argv[1]
c=cooler.Cooler(filename) ## replace with file names of interest
dic={chr:'chr'+str(chr) for chr in c.chromnames}
cooler.rename_chroms(c,dic)
print("Chromosome names updated successfully.")
