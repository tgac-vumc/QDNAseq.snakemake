bin <- commandArgs(TRUE)[1]
if (is.na(bin))
  bin <- 15
bin <- as.integer(bin)

genome <- commandArgs(TRUE)[2]
if (is.na(genome))
    genome <- "hg19"

library(QDNAseq)
library(Biobase)
library(R.cache)

setCacheRootPath(path="./.Rcache")

bins <- getBinAnnotations(bin, genome=genome)
QRC <- binReadCounts(bins, path='bam', cache=T)
saveRDS(QRC, paste(bin, "kbp-raw.rds", sep=""))
QRC.f <- applyFilters(QRC, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE)
QRC.f <- estimateCorrection(QRC.f, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE)
QRC.f <- applyFilters(QRC.f, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE, chromosome=NA)
QCN.fc <- correctBins(QRC.f)
QCN.fcn <- normalizeBins(QCN.fc)
QCN.fcns <- smoothOutlierBins(QCN.fcn)
saveRDS(QCN.fcns, paste(bin, "kbp.rds", sep=""))

QCN.fcnss <- segmentBins(QCN.fcns)
QCN.fcnssn <- normalizeSegmentedBins(QCN.fcnss)
saveRDS(QCN.fcnssn, paste(bin, "kbp-segmented.rds", sep=""))

QCN.fcnssnc <- callBins(QCN.fcnssn)
saveRDS(QCN.fcnssnc, paste(bin, "kbp-called.rds", sep=""))

exportBins(QCN.fcnssnc, paste(bin, "kbp-copynumbers.igv", sep=""), format="igv", type="copynumber")
exportBins(QCN.fcnssnc, paste(bin, "kbp-segments.igv", sep=""), format="igv", type="segments")
exportBins(QCN.fcnssnc, paste(bin, "kbp-calls.igv", sep=""), format="igv", type="calls")

