library(QDNAseq)
library(Biobase)

bam <- snakemake@input[[1]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
out <- snakemake@output[[1]]
genome <- snakemake@params[["genome"]]

#print(bam)
#print(bin)
#print(out)

bins <- getBinAnnotations(bin, genome=genome)
QRC <- binReadCounts(bins, bamfiles=bam, cache=F, chunkSize=10000000)

#sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples

#if (sum(duplicated(samples)) > 0) {
#        QRC <- poolRuns(QRC, samples=samples, force=TRUE)
#}

saveRDS(QRC, out)
#paste(out, "-", bin, "kbp-raw.rds", sep=""))
#QRC.f <- applyFilters(QRC, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE)
#QCN.fc <- correctBins(QRC.f)
#QCN.fcn <- normalizeBins(QCN.fc)
#QCN.fcns <- smoothOutlierBins(QCN.fcn)
#saveRDS(QCN.fcns, paste(bin, "kbp.rds", sep=""))

#QCN.fcnss <- segmentBins(QCN.fcns)
#QCN.fcnssn <- normalizeSegmentedBins(QCN.fcnss)
#saveRDS(QCN.fcnssn, paste(bin, "kbp-segmented.rds", sep=""))

#QCN.fcnssnc <- callBins(QCN.fcnssn)
#saveRDS(QCN.fcnssnc, paste(bin, "kbp-called.rds", sep=""))

#exportBins(QCN.fcnssnc, paste(bin, "kbp-copynumbers.igv", sep=""), format="igv", type="copynumber")
#exportBins(QCN.fcnssnc, paste(bin, "kbp-segments.igv", sep=""), format="igv", type="segments")
#exportBins(QCN.fcnssnc, paste(bin, "kbp-calls.igv", sep=""), format="igv", type="calls")
