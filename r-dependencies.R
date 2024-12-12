#!/usr/bin/env Rscript
##############################################################################

#Author: Erik Bosch
#date: Jan 20221

#This script is for installing remote packages to be used in the R-environment for snakemake
##############################################################################

library(devtools)
options(unzip = 'internal')

if (!require("BiocManager", quietly = TRUE))
	    install.packages("BiocManager",repos="http://cran.us.r-project.org")
BiocManager::install("DNAcopy")

# install WECCA from Github
devtools::install_github("tgac-vumc/WECCA", force=TRUE)

# install GGHtest from source
fn = 'http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtest_1.1.tar.gz'
if (!file.exists('CGHtest_1.1.tar.gz')){download.file(fn, destfile="CGHtest_1.1.tar.gz")} else {untar("CGHtest_1.1.tar.gz",list=TRUE)} # download or show content ofalready downloaded

install.packages("CGHtest_1.1.tar.gz", source=NULL, depenencies=TRUE)
