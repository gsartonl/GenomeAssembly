#!/bin/bash

#SBATCH --account pengel_spirit
#SBATCH --job-name qc_trimmomatic
#SBATCH --time 20:00
#SBATCH --mem 4G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/trimming.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/trimming.out

### Loading modules ###
module load gcc/9.3.0 # module required to run trimmomatic
module load trimmomatic/0.39 # trimmomatic module

### Variables ###
workingDir=/users/gsartonl/spirit/gsartonl/assembly/
inDir=${workingDir}/00_rawReads
outDir=${workingDir}/01_trimming
yourfileID=10_S10

### commands ###
mkdir -p ${outDir}

trimmomatic PE -threads 8 ${inDir}/${yourfileID}_L001_R1_001.fastq.gz ${inDir}/${yourfileID}_L001_R2_001.fastq.gz \
${outDir}/${yourfileID}_R1_paired.fastq.gz ${outDir}/${yourfileID}_R1_unpaired.fastq.gz \
${outDir}/${yourfileID}_R2_paired.fastq.gz ${outDir}/${yourfileID}_R2_unpaired.fastq.gz \
ILLUMINACLIP:NexteraPE-PE.fa:3:25:6 LEADING:9 TRAILING:9 \
SLIDINGWINDOW:4:15 MINLEN:60 

### END ###
