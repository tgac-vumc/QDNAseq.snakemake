library(QDNAseq)
library(Biobase)

if(!require(BSgenome.Hsapiens.UCSC.hg38)) {
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
}

bigWig <- snakemake@input[[1]]
binSize <- as.integer(snakemake@params[["binSize"]])

print(bigWig)
print(binSize)

# Create bin annotations
bins <- createBins(BSgenome.Hsapiens.UCSC.hg38, binSize)

# Calculate mappability bin average
bins$mappability <- calculateMappability(bins, bigWigFile=bigWig, bigWigAverageOverBed='bigWigAverageOverBed')
saveRDS(bins, paste("rds/mappabilityBins-", binSize, "k", ".rds", sep=""))
