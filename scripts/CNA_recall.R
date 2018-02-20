#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis segmentBins
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(CGHcall))
suppressMessages(library(CGHtest))
suppressMessages(library(denstrip))

source('scripts/CGHcallPlus.R')
source('scripts/reCall5levelCNAs.R')
source('scripts/plotQDNAseq.R')

called <- snakemake@input[["called"]]
recalled <- snakemake@output[["recalled"]]
profiles <- snakemake@params[["profiles"]]
copynumbers<-snakemake@output[["copynumbers"]]
segments<-snakemake@output[["segments"]]
calls<-snakemake@output[["calls"]]
copynumbersbed<-snakemake@params[["copynumbersbed"]]
segmentsbed<-snakemake@params[["segmentsbed"]]
failed <- snakemake@params[["failed"]]
bedfolder <- snakemake@params[["bedfolder"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=TRUE, split=FALSE)

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
exportBins(reCalledRCs, file=copynumbersbed, format="bed", logTransform=TRUE, type="copynumber")
exportBins(reCalledRCs, file=segmentsbed, format="bed", type="segments")

#create output for failed samples - for snakemake compatibility.
failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){for(file in failed_samples[,1]){file.create(paste(profiles, file,".png",sep=""))
file.create(paste(bedfolder, file,"-copynumbers.bed",sep=""))
file.create(paste(bedfolder, file,"-segments.bed",sep=""))
}}
