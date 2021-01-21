#!/usr/bin/env Rscript
##############################################################################

#Author: Erik Bosch
#date: Jan 20221

#This script is for installing remote packages to be used in the R-environment for snakemake
# snakemake --use-conda --cores=1 -U create_r_environment --rerun-incomplete
##############################################################################

library(devtools)
options(unzip = 'internal')
#devtools::install_version("denstrip", version = '1.5.4')
#devtools::install_version("flexmix", version = '2.3-17')
#devtools::install_version("flexmix"')
#install.packages("flexmix")

devtools::install_github("tgac-vumc/WECCA")

#install.packages(file.path(getwd(),"CGHtest_1.1.tar.gz"))
install.packages("CGHtest_1.1.tar.gz", source=NULL, depenencies=TRUE)

