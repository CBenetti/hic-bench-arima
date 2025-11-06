library(HiCExperiment)
library(HiContacts)
library(BiocParallel)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(rtracklayer)


pc1 <- import("results/compartments.by_sample.homer.res_100kb/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/CELL-HEMATO-NK-INTER_S3/compartments.scores.bedGraph")
seqlevels(pc1) <- gsub("^chr", "", seqlevels(pc1))
seqnames(pc1)<-gsub("^chr", "", as.character(seqnames(pc1)))
hic_file<-"../tracks/results/tracks.by_sample.juicer/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/CELL-HEMATO-NK-INTER_S3/filtered.hic"
hic_exp<-import(hic_file,format="hic",resolution=100000)
common <- intersect(seqlevels(hic_exp), seqlevels(pc1))
pc1 <- keepSeqlevels(pc1, common, pruning.mode = "coarse")
colnames(mcols(pc1)) <- "PC1"
topologicalFeatures(hic_exp, "compartments")<- as(pc1, "GRanges")
metadata(hic_exp)$eigens <- as(pc1, "GRanges")
mcols(metadata(hic_exp)$eigens)$eigen<- mcols(pc1)$PC1
mcols(topologicalFeatures(hic_exp,"compartments"))$compartment <- ifelse(mcols(topologicalFeatures(hic_exp,"compartments"))$PC1>0,"A","B")
pdf("results/compartments.by_sample.homer.res_100kb/filter.by_sample.mapq_30_mindist0/align.by_sample.hicpro/CELL-HEMATO-NK-INTER_S3/saddle_plot.pdf")
plotSaddle(hic_exp, nbins = 25, BPPARAM = SerialParam(progressbar = FALSE))
dev.off()
