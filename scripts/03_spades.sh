#!/bin/bash

# script for SPAdes de novo assembly #
## @ gsartonl ##
## 07.01.2022 ##

### Variables ###
FWD_trimmed=$1
REV_trimmed=$2
outDir=$3
threads=$4

### commands ###
spades.py -1 ${FWD_trimmed} -2 ${REV_trimmed} --careful --threads ${threads} -o ${outDir}

### END ###
