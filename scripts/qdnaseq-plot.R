library(QDNAseq)
library(Biobase)
library(matrixStats)

input <- commandArgs(TRUE)[1]
genomes <- commandArgs(TRUE)[2]
chrs <- commandArgs(TRUE)[3]
if (is.na(input))
  stop('No input file provided.\n')
if (is.na(genomes))
  genomes <- 'yes'
if (is.na(chrs))
  chrs <- 'no'

output <- sub('\\.rds$', '', input)
output <- sub('-normalized', '', output)
output <- sub('-segmented', '', output)
output <- sub('-called', '', output)

bin <- output

dat <- readRDS(input)
if ('calls' %in% assayDataElementNames(dat)) {
  output <- paste0(output, '-called')
  png(paste0(output, ".png"), width=297, height=210, units='mm', res=150)
#  pnga4(paste0(output, '.png'))
  frequencyPlot(dat)
  dev.off()
} else if ('segmented' %in% assayDataElementNames(dat)) {
  output <- paste0(output, '-segmented')
} else {
  output <- paste0(output, '-normalized')
}

plot.profiles <- function(cgh, directory, byChr=FALSE) {
  tmp <- sampleNames(cgh)
  if (!file.exists(directory))
    dir.create(directory)
  if ('filter' %in% colnames(fData(cgh))) {
    chrs <- unique(chromosomes(cgh)[fData(cgh)$filter])
  } else {
    chrs <- unique(chromosomes(cgh))
  }
  for (i in 1:length(sampleNames(cgh))) {
    if (byChr) {
      for (chr in chrs) {
        png(file.path(directory, paste(tmp[i], '-chr', chr, '.png', sep='')), width=297, height=210, units='mm', res=150)
	print(sum(chromosomes(cgh) == chr))
	Sys.sleep(1)
        plot(cgh[chromosomes(cgh) == chr,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
        dev.off()
      }
    } else {
      png(file.path(directory, paste(tmp[i], '.png', sep='')), width=297, height=210, units='mm', res=150)
      plot(cgh[,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
      dev.off()
    }
  }
}

if (file.exists('stats')) {
  mads <- apply(assayDataElement(dat, 'copynumber'), 2, madDiff, na.rm=TRUE)
  for (s in sampleNames(dat))
    cat(sprintf("%.3f", mads[s]), file=file.path('stats', paste(s, '.mad.', bin, sep="")))
}


if (genomes=='yes')
  plot.profiles(dat, output)
if (chrs=='yes')
  plot.profiles(dat, paste(output, '-bychr', sep=""), byChr=TRUE)



# EOF
