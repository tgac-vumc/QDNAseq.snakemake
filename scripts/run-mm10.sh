#!/bin/bash

# Create fastq directory and make all data sets available underneath this dir
# Pipeline will recursively find all .fastq or fastq.gz files and process them
 
# Pipeline dir
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export ref="/ccagc/data/ref/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.5.x/genome.fa"
export genome="mm10"

$DIR/bwa-align-all.sh
