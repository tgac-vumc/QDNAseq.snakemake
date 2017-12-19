#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis to call bins
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(CGHcall))
suppressMessages(library(CGHtest))

source('scripts/CGHcallPlus.R')
source('scripts/plotQDNAseq.R')

segmented <- snakemake@input[["segmented"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
called <- snakemake@output[["called"]]
profiles <- snakemake@params[["profiles"]]
freqplots <- snakemake@output[["freqplot"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)

##############################################################################################################
# Call data & frequency plot calls
##############################################################################################################

# load data
QCN.fcnsdsn <- readRDS(segmented)

# perform functions
# Call withouth correction for cellularity
QCN.fcnsdsnc <- callBins(QCN.fcnsdsn, nclass=5)
saveRDS(QCN.fcnsdsnc, called)

# Plot called profiles
plotQDNAseq(QCN.fcnsdsnc, profiles)

# frequency plot calls
png(freqplots, res=300, width=14, height=7, unit='in')
frequencyPlot(QCN.fcnsdsnc)
dev.off()
