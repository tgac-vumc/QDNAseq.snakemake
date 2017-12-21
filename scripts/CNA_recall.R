#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis segmentBins
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
library(QDNAseq)
library(Biobase)
library(CGHcall)
library(CGHtest)
library(denstrip)

source('scripts/CGHcallPlus.R')
source('scripts/reCall5levelCNAs.R')
source('scripts/plotQDNAseq.R')

called <- snakemake@input[["called"]]
recalled <- snakemake@output[["recalled"]]
profiles <- snakemake@params[["profiles"]]
copynumbers<-snakemake@output[["copynumbers"]]
segments<-snakemake@output[["segments"]]
calls<-snakemake@output[["calls"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)

##############################################################################################################
# RECALL data
##############################################################################################################
# load data
QCN.fcnsdsnc <- readRDS(called)

# Re-call: copyNumbersCalled in 5 level CNAs
reCalledRCs <- reCall5levelCNAs(QCN.fcnsdsnc)
saveRDS(reCalledRCs, recalled)
# Plot re-called profiles
plotQDNAseq(reCalledRCs, profiles)

##############################################################################################################
# Create IGV objects from readcounts
##############################################################################################################

exportBins(reCalledRCs, copynumbers, format="igv", type="copynumber")
exportBins(reCalledRCs, segments, format="igv", type="segments")
exportBins(reCalledRCs, calls, format="igv", logTransform=FALSE, type="calls")
