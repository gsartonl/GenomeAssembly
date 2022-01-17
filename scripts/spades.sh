#!/bin/bash

#SBATCH --account pengel_spirit
#SBATCH --job-name spades
#SBATCH --time 45:00
#SBATCH --mem 8G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/assembly.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/assembly.out

### Loading modules ###
module load gcc/9.3.0 # module required to run spades
module load python    # module required to run spades
module load fastqc/0.11.9 # spades module

### Variables ###
workingDir=/users/gsartonl/spirit/gsartonl/assembly/
inDir=${workingDir}/01_trimming
outDir=${workingDir}/02_assembly

yourfileID=10_S10

### commands ###

mkdir -p ${outDir} #create a directory in which the fastqc output will be placed

spades.py -1 ${inDir}/${yourfileID}_R1_paired.fastq.gz -2 ${inDir}/${yourfileID}_R2_paired.fastq.gz --careful -o ${outDir} 

### END ###