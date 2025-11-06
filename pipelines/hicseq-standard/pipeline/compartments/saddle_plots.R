library(HiCExperiment)
library(HiContacts)
library(BiocParallel)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(rtracklayer)

# --- Read sample name from command line ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript saddle_plots.R <SAMPLE_NAME>")
sample <- args[1]

# --- Build dynamic paths ---
pc1_path <- sprintf("results/compartments.by_sample.homer.res_100kb/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/%s/compartments.scores.bedGraph", sample)
hic_file <- sprintf("../tracks/results/tracks.by_sample.juicer/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/%s/filtered.hic", sample)
pdf_out <- sprintf("results/compartments.by_sample.homer.res_100kb/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/%s/saddle_plot.pdf", sample)

# --- Run analysis ---
pc1 <- import(pc1_path)
seqlevels(pc1) <- gsub("^chr", "", seqlevels(pc1))
seqnames(pc1) <- gsub("^chr", "", as.character(seqnames(pc1)))
hic_exp <- import(hic_file, format="hic", resolution=100000)

common <- intersect(seqlevels(hic_exp), seqlevels(pc1))
pc1 <- keepSeqlevels(pc1, common, pruning.mode="coarse")
colnames(mcols(pc1)) <- "PC1"
topologicalFeatures(hic_exp, "compartments") <- as(pc1, "GRanges")
metadata(hic_exp)$eigens <- as(pc1, "GRanges")
mcols(metadata(hic_exp)$eigens)$eigen <- mcols(pc1)$PC1
mcols(topologicalFeatures(hic_exp,"compartments"))$compartment <- ifelse(
  mcols(topologicalFeatures(hic_exp,"compartments"))$PC1 > 0, "A", "B"
)

pdf(pdf_out)
plotSaddle(hic_exp, nbins = 25, BPPARAM = SerialParam(progressbar = FALSE))
dev.off()
