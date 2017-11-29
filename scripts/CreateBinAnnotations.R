library(QDNAseq)
library(Biobase)

binSize <- snakemake@input[[1]]

# Load mappability bins
bins <- readRDS(paste("rds/mappabilityBins-",binSize,"k",".rds",sep="")) -> bins

# Calculate blacklist 
# bins$blacklist <- calculateBlacklist(bins,bedFiles=c("hg38.blacklist.bed"))

# Calculate residuals
ctrl <- binReadCounts(bins,path="bam")
saveRDS(ctrl, paste("rds/1kg.readCounts.", binSize, "k", ".rds", sep=""))

ctrl <- applyFilters(ctrl, residual=FALSE, blacklist=FALSE,mappability=FALSE, bases=FALSE)

bins$residual <- iterateResiduals(ctrl)

bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0

bins <- AnnotatedDataFrame(bins,
varMetadata=data.frame(labelDescription=c(
'Chromosome name',
'Base pair start position',
'Base pair end position',
'Percentage of non-N nucleotides (of full bin size)',
'Percentage of C and G nucleotides (of non-N nucleotides)',
'Average mappability of 50mers with a maximum of 2 mismatches',
'Percent overlap with ENCODE blacklisted regions',
'Median loess residual from 1000 Genomes (50mers)',
'Whether the bin should be used in subsequent analysis steps'),
row.names=colnames(bins)))

attr(bins, "QDNAseq") <- list(
author="Daoud Sie",
date=Sys.time(),
organism="Hsapiens",
build="hg38",
version=packageVersion("QDNAseq"),
md5=digest::digest(bins@data),
sessionInfo=sessionInfo())

saveRDS(bins,paste("rds/FinalBins-",binSize,"k",".rds",sep=""))
