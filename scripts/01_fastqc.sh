#!/bin/bash

####-------------------------------------------------####
# Reads trimming for genome assembly with SNAKEMAKE #
##           author : gsartonl          ##
##           date : 06/01/2022           ##
####-------------------------------------------------####

### Variables ###

FWD=$1;
REV=$2;
outDir=$3;

### Output directory ###
mkdir -p ${outDir}

### commands ###

fastqc --outdir ${outDir} ${FWD} ${REV}

### END ###
