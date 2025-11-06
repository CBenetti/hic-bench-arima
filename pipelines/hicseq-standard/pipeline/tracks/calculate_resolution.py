import sys
import cooler
import pandas as pd
import numpy as np
import os

# ---- 1. Get sample name from command-line argument ----
if len(sys.argv) < 3:
    print("Usage: python compute_coverage.py SAMPLE_NAME BRANCH")
    sys.exit(1)

sample = sys.argv[1]
base_dir = sys.argv[2] #e.g. "results/tracks.by_sample.juicer_cool.multires/filter.by_sample.mapq_20_mindist0/align.by_sample.hicpro"

# ---- 2. Base paths ----
#base_dir = "results/tracks.by_sample.juicer_cool.multires/filter.by_sample.mapq_20_mindist0/align.by_sample.hicpro"
sample_dir = os.path.join(base_dir, sample)
mcool_path = os.path.join(sample_dir, "filtered.mcool")

# ---- 3. Loop over resolutions ----
for res in [40000, 100000]:
    print(f"\nProcessing {sample} at {res:,} bp resolution...")

    # Open cooler object
    uri = f"{mcool_path}::resolutions/{res}"
    c = cooler.Cooler(uri)

    # Extract bins and pixels
    bins = c.bins()[:][['chrom', 'start', 'end']]
    pixels = c.pixels()[:]

    # Compute per-bin coverage
    cov = (
        pd.concat([
            pixels.groupby('bin1_id')['count'].sum(),
            pixels.groupby('bin2_id')['count'].sum()
        ])
        .groupby(level=0)
        .sum()
    )

    cov = cov.reindex(range(len(bins)), fill_value=0)  # ensures all bins are covered
    bins['coverage'] = cov.fillna(0)

    # ---- 4. Compute 20th percentile ----
    perc20 = np.nanpercentile(bins['coverage'].dropna(), 20)
    perc20_rounded = round(perc20, 2)

    # ---- 5. Save output with percentile in filename ----
    out_file = os.path.join(
        sample_dir,
        f"coverage_{res//1000}kb_p20-{perc20_rounded}.tsv"
    )
    bins.to_csv(out_file, sep="\t", index=False)
