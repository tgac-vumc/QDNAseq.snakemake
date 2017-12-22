#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis
# Function to create a BED file of all regions per sample
## Date: 8 August 2017, December 2017
# Author: Matias Mendeville
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(CGHregions))

source("scripts/addCytobands.R")
source("scripts/CGHcallPlus.R")

RegionsCGH<-snakemake@input[["RegionsCGH"]]
max.focal.size.mb<-snakemake@params[["max_focal_size_mb"]]
min.freq.focal<-snakemake@params[["min_freq_focal"]]
allRegions<-snakemake@output[["allRegions"]]
allFocalRegions<-snakemake@output[["allFocalRegions"]]
cytoband_data<-snakemake@params[["cytobands"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##############################################################################################################
# Create CGHregion BED file + Annotate
##############################################################################################################

# load data
QCN.Regions <- readRDS(RegionsCGH)

# create BED file from CGHregions

makeCGHregionsTable <- function(CGHregions, max.focal.size.mb=3, min.freq.focal=25, output_allRegions, output_allFocalRegions,  cytoband_data){

	# all regions in BED format
	regionsBED <- fData(CGHregions)[,1:3]

	# region sizes
	regionSize <- fData(CGHregions)[,3] - fData(CGHregions)[,2]
	regionsBED <- cbind(regionsBED, regionSize)

	# add cytoband info
	regionsBED <- addCytobands(regionsBED, cytoband_data)

	# frequencies of gains and losses
	calls <- regions(CGHregions)
	loss.freq <- rowMeans(calls < 0)
	gain.freq <- rowMeans(calls > 0)

	# add to 1 table
	regionData <- cbind(regionsBED, loss.freq, gain.freq)
	colnames(regionData) <- c('chromo', 'start', 'end', 'Region size (Mb)', 'Cytoband', 'loss (%)', 'gain (%)')

	# round integers:
	regionData[,4] <- round(regionData[,4]/1000000, digits=2)
	regionData[,c(6,7)] <- round(regionData[,c(6,7)]*100, digits=2)

	# save
	write.table(regionData, output_allRegions, sep='\t', row.names=F, col.names=T, quote=F)

	# make separate table for focals only
	focal.ix <- which(regionData[,4] <= max.focal.size.mb)
	focalRegions <- regionData[focal.ix, ]

	# and in more than 25% of samples
	loss.cutoff.ix <- which(focalRegions[,6] > min.freq.focal)
	gain.cutoff.ix <- which(focalRegions[,7] > min.freq.focal)
	cutoff.ix <- sort( c(loss.cutoff.ix, gain.cutoff.ix ))

	focalRegions <- focalRegions[cutoff.ix, ]
	if(nrow(focalRegions) != 0 ){rownames(focalRegions) <- c(1:nrow(focalRegions))}

	write.table(focalRegions, output_allFocalRegions, sep='\t', row.names=F, col.names=F, quote=F)
}

makeCGHregionsTable(QCN.Regions, max.focal.size.mb, min.freq.focal, allRegions, allFocalRegions,  cytoband_data)

# annotate BED files
#system("bash $DIR/annotateFocalCNAbed.sh")
#system("bash /ccagc/home/matias/tools/MM_QDNAseq/annotateRegionsFocalCNAbed.sh")
