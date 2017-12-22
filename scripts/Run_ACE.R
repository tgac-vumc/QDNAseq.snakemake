#!/usr/bin/env Rscript
##############################################################################

#Author: Tjitske logs
#date: Dec 2017

#This script is a small wrapper around ACE to work in the snakemake pipeline
##############################################################################

suppressMessages(library(QDNAseq))
source('scripts/ACE.R')

ploidies<-as.integer(snakemake@wildcards[["ploidy"]])
binsizes<-as.integer(snakemake@wildcards[["binSize"]])
inputfile <-snakemake@input[["segmented"]]
outputdir<-snakemake@params[["outputdir"]]
log<-snakemake@log[[1]]

imagetype <- snakemake@config[["ACE"]][["imagetype"]]
method<-snakemake@config[["ACE"]][["method"]]
penalty<-as.numeric(snakemake@config[["ACE"]][["penalty"]])
cap<-as.integer(snakemake@config[["ACE"]][["cap"]])
trncname<- snakemake@config[["ACE"]][["trncname"]]
printsummaries<- snakemake@config[["ACE"]][["printsummaries"]]

copyNumbersSegmented <- readRDS(inputfile)

parameters <- data.frame(options = c("ploidies","imagetype","method","penalty","cap","trncname","printsummaries"),
                         values = c(paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,trncname,printsummaries))

write.table(parameters, file=log, quote = FALSE, sep = "\t", na = "", row.names = FALSE) 

ploidyplotloop(copyNumbersSegmented ,outputdir , ploidies,imagetype,method,penalty,cap,trncname,printsummaries)
