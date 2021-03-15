msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(CGHtest))
} else{
library(CGHtest)
}

source("scripts/CGHtest_functions.R", echo = FALSE)

##########
# input file paths
CGHregions_RDS_file_path <-"./output-CGH/100kbp/data/100kbp-RegionsCGH.rds"
CGHregions_RDS_file_path <- snakemake@input[["RegionsCGH"]]

#output file paths
freqPlotCompare_output <-'./output-CGH/CGHtest/freqPlotCompare.png' 
CGHtest_output <-'./output-CGH/CGHtest/CGHtest.Rdata'
plotPFDR_output <-'./output-CGH/CGHtest/plotPFDR.png'
combined_output <-'./output-CGH/CGHtest/combined.png'

freqPlotCompare_output <- snakemake@output[["freqPlotCompare"]]
CGHtest_output <- snakemake@output[["CGHtest"]]
plotPFDR_output <- snakemake@output[["plotPFDR"]]
combined_output <- snakemake@output[["combined"]]

# params path
outputdir <- snakemake@params[["outputdir"]]

print(CGHregions_RDS_file_path)
print(freqPlotCompare_output)
print(plotPFDR_output)
print(combined_output)



log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

#if(!dir.exists('./output-CGH/CGHtest')){dir.create('./output-CGH/CGHtest')}
if(!dir.exists(outputdir)){dir.create(outputdir)}

# read CGHregions RDS
CGHregions_RDS_file<- readRDS(CGHregions_RDS_file_path)

# define groups to compare
group1_index <- 1:9
group2_index <- 10:18
  
# subset data
group1 <- CGHregions_RDS_file[,group1_index]
group2 <- CGHregions_RDS_file[,group2_index]

# Compare 2 cohorts
teststat <- "Chi-square"
compareCNAs2cohorts(x = CGHregions_RDS_file,
                    data1 = group1,
                    data2 = group2,
                    cohort1_ID = "group1",
                    cohort2_ID = "group2",
                    teststat,
                    output_freqPlotCompare = freqPlotCompare_output,
                    output_CGHtest = CGHtest_output,
                    output_plotPFDR = plotPFDR_output,
                    output_combined = combined_output)
