#!/bin/bash

##  Used to generate test dataset of raw fastq files (00-raw-fastq)
##  Takes input fastq file and splits it into two (to simulate the DNAseq core raw reads format)
##  Usage: split_fastq_SE.sh <output_dir>
##  NOTE: For it to work, run this script in the directory with the de-multiplexed ".gz" fastq read files
##  Currently for single-end reads only
##  abhasi / July 25, 2016

ARGS=1
if [ $# -ne "$ARGS" ]
then
  echo "Usage: `basename $0` output_dir";
  exit $E_BADARGS;
fi

for sample in Sample*
  do 
    mkdir -p $1/${sample/_R1*/};
    length="$(zcat $sample | wc -l)"; 
    half="$((${length/ Sample_*/} / 2))";
    
    echo "${sample/_R1*/} (length: ${length/ Sample_*/}) was split into 2 files, each with ${half} lines";
    zcat $sample | head -n $half > $1/${sample/_R1*/}/${sample/_R1*/_S16_L007_001_R1.fastq};
    zcat $sample | tail -n $half > $1/${sample/_R1*/}/${sample/_R1*/_S16_L007_002_R1.fastq};
    gzip $1/${sample/_R1*/}/* ;
    
   done
