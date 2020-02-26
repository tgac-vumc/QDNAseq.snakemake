#First part is to call with altered cellularity using chgCall.
#The easiest is to do everything from the QDNAseq-snakemake folder and with the conda environment because than you are sure that the same programversions are used.

#from the terminal

#cd /net/nfs/PAT/analysis/MPS-343/Snakemake/QDNAseq.snakemake

#source activate QDNAseq-snakemake

#Open R

#R

#open required libraries.

library(QDNAseq)
library(Biobase)
library(CGHcall)
library(CGHtest)

source('scripts/CGHcallPlus.R')
source('scripts/plotQDNAseq.R')

segmented <- "../100kbp/data/100kbp-segmented.rds"
bin<-100
called <- "../100kbp/data/100kbp-called-cellularity.rds"
freqplots<- "../100kbp/profiles/freqPlot/freqPlot_100kbp_cellularity.png"
profiles <- "../100kbp/profiles/called_cellularity/"
fitpickertable<- "../100kbp/ACE/2N/fitpicker_2N.tsv"
ACE_output<-"../100kbp/ACE/2N/postanalysis/"

##############################################################################################################
# Call data & frequency plot calls
##############################################################################################################

# load data
QCN.fcnsdsn <- readRDS(segmented)
fitpicker<- read.delim(fitpickertable, stringsAsFactors=FALSE)

#perform functions
# Call with correction for cellularity
QCN.fcnsdsnc <- callBins(QCN.fcnsdsn, nclass=5, cellularity=fitpicker$likely_fit)
saveRDS(QCN.fcnsdsnc, called)

dir.create(profiles)
plotQDNAseq(QCN.fcnsdsnc, profiles)

# frequency plot calls
png(freqplots, res=300, width=14, height=7, unit='in')
frequencyPlot(QCN.fcnsdsnc)
dev.off()
