#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
genome=${genome:-"hg19"}

module load QDNAseq/1.12
module load FastQC

if [ ! -d fastq ]
then
  echo ERROR: fastq directory not found.
  exit
fi

if [ ! -d qc-fastq ]
then
  mkdir qc-fastq
fi

if [ ! -d bam ]
then
  mkdir bam
fi

if [ ! -d stats ]
then
  mkdir stats
fi

if [ ! -d qc-bam ]
then
  mkdir qc-bam
fi

cd bam

files=`find -L ../fastq/ -type f | grep -e "fastq$" -e "fastq.gz$" -e "fq$" -e "fq.gz$"`

for file in $files
do
  f=`basename $file .gz`
  f=`basename $f .fastq`
  f=`basename $f .fq`
  d=`dirname $file`
  echo `date` " - $f Alignment started"
  if [ -f $f.bam ]
  then
    continue
  fi
  $DIR/bwa-align.sh $d/$f
  echo `date` " - $f Alignment done"
done

cd ..

$DIR/putLightBox.sh ./

for bin in 15 30 100 1000
do
  if [ ! -e ${bin}kbp-raw.rds ]
  then
    Rscript $DIR/QDNAseq.R $bin $genome
  fi
  Rscript $DIR/qdnaseq-plot.R ${bin}kbp.rds yes no
  Rscript $DIR/qdnaseq-plot.R ${bin}kbp-segmented.rds yes no
  for i in `find ./ -maxdepth 1 -type "d" -iname "${bin}kbp*"`
  do
    cd $i	
    $DIR/createLightBox.sh > index.html
    cd ..
  done
done

$DIR/lane-summary.sh > summary.html

fqs=`find -L ./fastq/ -type f | grep -e "fastq$" -e "fastq.gz$" -e "fq$" -e "fq.gz$"`
fastqc $fqs -o qc-fastq && rm qc-fastq/*_fastqc.zip
fastqc bam/*.bam -f bam_mapped -o qc-bam && rm qc-bam/*_fastqc.zip

# EOF
