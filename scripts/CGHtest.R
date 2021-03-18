###########
# derived by Erik Bosch from the script of Jurriaan Janssen, original from Matias Mendeville 16-03-2021
###########

msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(CGHtest))
suppressMessages(library(readxl))
} else{
library(CGHtest)
library(readxl)
}

source("scripts/CGHtest_functions.R", echo = FALSE)

##########
# input file paths
CGHregions_RDS_file_path = '../../../output-CGH/100kbp/data/100kbp-RegionsCGH.rds'
CGHregions_RDS_file_path <- snakemake@input[["RegionsCGH"]]

freqPlotCompare_output <- snakemake@output[["freqPlotCompare"]]
CGHtest_output <- snakemake@output[["CGHtest"]]
plotPFDR_output <- snakemake@output[["plotPFDR"]]
combined_output <- snakemake@output[["combined"]]

# params path
outputdir <- snakemake@params[["outputdir"]]
clinicaldataPath <- 'clinical.xlsx'
clinicaldataPath <- snakemake@params[["clinicaldataPath"]]

clinicaldata_col_snames <- 1
clinicaldata_col_snames <- snakemake@params[["columnSampleNames"]]

set_clinicaldata_class_samples <- c('short','long')
set_clinicaldata_class_samples <- snakemake@params[["ClassSamples"]]
clinicaldata_col_class_samples <- 3
clinicaldata_col_class_samples <- snakemake@params[["columnClassSamples"]]


log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##################

if(!dir.exists(outputdir)){dir.create(outputdir)}

# read CGHregions RDS
CGHregions_RDS_file <- readRDS(CGHregions_RDS_file_path)

# get sample names from CGHregions
snames <- sampleNames(CGHregions_RDS_file)

# read clinical data
clinicaldata <- as.data.frame(read_excel(clinicaldataPath))
clinicaldata_snames  <- clinicaldata[,clinicaldata_col_snames]

if (!isTRUE(all.equal(sort(snames),sort(clinicaldata_snames)))){
    print(data.frame('Clinical Sample Names' = clinicaldata_snames, 'CGHregions Sample Names' = snames))
    stop('Sample Names in clinical data and CGHregions RDS file do not match!')
}
# match sample names from CGHregions with sample names in clinical data
snames_i  <- match(clinicaldata_snames, snames)

# redefine clinical data matching the CGHregions names
clinicaldata <- clinicaldata[snames_i,]

# define clinical data classes and match with set classes
clinicaldata_class_samples_all <- clinicaldata[,clinicaldata_col_class_samples]
clinicaldata_class_samples_unique <- unique(clinicaldata_class_samples_all)

tempa <- clinicaldata_class_samples_unique
tempb <- set_clinicaldata_class_samples
if (!isTRUE(all.equal(sort(tempa),sort(tempb)))){
    print(data.frame('Classes in clinical data' = tempa, 'Classes defined by user in config' = tempb))
    stop('Classes in the clinical data do not match the classes defined for the CGHtest-setting in the configuration-file!')
}

# define groups to compare
group1_index <- which(clinicaldata_class_samples_all == set_clinicaldata_class_samples[1]) # c(4,5,8,11,14,17) for clinical.xlsx
group2_index <- which(clinicaldata_class_samples_all == set_clinicaldata_class_samples[2]) # c(1,2,3,6,7,9,10,12,13,15,16,18) for clinical.xlsx
  
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
