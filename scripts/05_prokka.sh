#!/bin/bash


echo "#####################################################"
echo $id "Starting time: $(date -u)"
echo "#####################################################"


# ### Create directories and define variables
workdir=$1
assembly=$2
outDir=$3;
genomeID=$4;

### command

prokka --outdir ${outDir} --locustag ${genomeID} --prefix ${genomeID} ${assembly} 
