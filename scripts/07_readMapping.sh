#!/bin/bash

### Script for ReadMapping ##

assembly=${1}/contigs.fasta
FWD=${2};
REV=${3};
genomeID=${4}
sam=$5;
bam=$6;
sorted=$7
depth=$8

bowtie2-build -f ${assembly} ${genomeID} --threads 16

bowtie2 -x ${genomeID} -1 ${FWD} -2 ${REV} --sensitive > ${sam}

samtools faidx ${assembly}

samtools view -S -b ${sam} > ${bam}

samtools sort ${bam} > ${sorted}

samtools depth -a  ${sorted} > ${depth}
