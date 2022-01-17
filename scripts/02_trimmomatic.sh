#!/bin/bash



### Variables ###

echo $1 $2 $3 $4

FWD=$1;
REV=$2;
outDir=$3;
genomeID=$4;
adaptors=$5;




### commands ###
mkdir -p ${outDir}/unpaired

trimmomatic PE -threads 8 ${FWD} ${REV} \
${outDir}/${genomeID}_R1_paired.fastq.gz ${outDir}/unpaired/${genomeID}_R1_unpaired.fastq.gz \
${outDir}/${genomeID}_R2_paired.fastq.gz ${outDir}/unpaired/${genomeID}_R2_unpaired.fastq.gz \
ILLUMINACLIP:${adaptors}:3:25:6 LEADING:9 TRAILING:9 \
SLIDINGWINDOW:4:15 MINLEN:60

### END ###
