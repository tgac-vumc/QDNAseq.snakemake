#!/usr/bin/env Rscript
##############################################################################################################
# script for calling based on cellularity estimated with ACE
# date: October 2018
# Tjitske Los
##############################################################################################################
#open required libraries.
library(QDNAseq)
library(Biobase)
library(CGHcall)
library(CGHtest)
#suppressMessages(library(QDNAseq))
#suppressMessages(library(Biobase))
#suppressMessages(library(CGHcall))
#suppressMessages(library(CGHtest))
#suppressMessages(library(denstrip))

source('scripts/CGHcallPlus.R')
source('scripts/plotQDNAseq.R')

called <- snakemake@input[["called"]]
fitpickertable<- snakemake@input[["fitpicker"]]

recalled <- snakemake@output[["recalled"]]
#copynumbers<-snakemake@output[["copynumbers"]]#moved to QDNAseq-segment
#segments<-snakemake@output[["segments"]]  #moved to QDNAseq-segment
calls<-snakemake@output[["calls"]]
profiles <- snakemake@params[["profiles"]]
#copynumbersbed<-snakemake@params[["copynumbersbed"]]  #moved to QDNAseq-segment
#segmentsbed<-snakemake@params[["segmentsbed"]] #moved to QDNAseq-segment
failed <- snakemake@params[["failed"]]
#bedfolder <- snakemake@params[["bedfolder"]]
freqplot<- snakemake@output[["freqplot"]]

minimum_cellularity <- snakemake@params[["minimum_cellularity"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=TRUE, split=FALSE)

##############################################################################################################
# Call data & frequency plot calls
##############################################################################################################
#First part is to call with altered cellularity using chgCall.
# load data
QCN.fcnsdsn <- readRDS(called)
fitpicker<- read.delim(fitpickertable, stringsAsFactors=FALSE)
fitpicker$likely_fit[fitpicker$likely_fit < minimum_cellularity] <- minimum_cellularity

#perform functions
# Call with correction for cellularity
reCalled <- callBins(QCN.fcnsdsn, nclass=5, cellularity=fitpicker$likely_fit)
saveRDS(reCalled, recalled)

##############################################################################################################
# Create profiles and frequency plot
##############################################################################################################
fitpicker$likely_fit<-paste0("cellularity: ",fitpicker$likely_fit )

dir.create(profiles)
plotQDNAseq(reCalled, profiles , subtitle=fitpicker$likely_fit)

# frequency plot calls
png(freqplot, res=300, width=14, height=7, unit='in')
frequencyPlot(reCalled)
dev.off()

##############################################################################################################
# Create IGV objects from readcounts
##############################################################################################################

#exportBins(reCalled, copynumbers, format="igv", type="copynumber") moved to CNA segments
#exportBins(reCalled, segments, format="igv", type="segments") moved to CNA segments
exportBins(reCalled, calls, format="igv", logTransform=FALSE, type="calls")
#exportBins(reCalled, file=copynumbersbed, format="bed", logTransform=TRUE, type="copynumber")
#exportBins(reCalled, file=segmentsbed, format="bed", type="segments")

#create output for failed samples - for snakemake compatibility.
failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){for(file in failed_samples[,1]){file.create(paste(profiles, file,".png",sep=""))
#file.create(paste(bedfolder, file,"-copynumbers.bed",sep=""))
#file.create(paste(bedfolder, file,"-segments.bed",sep=""))
}}
