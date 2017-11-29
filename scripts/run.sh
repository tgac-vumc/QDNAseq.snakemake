#!/bin/bash

# Create fastq directory and make all data sets available underneath this dir
# Pipeline will recursively find all .fastq or fastq.gz files and process them
 
# Pipeline dir
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

$DIR/bwa-align-all.sh
