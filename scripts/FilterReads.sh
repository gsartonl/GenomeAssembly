#!/bin/bash


#SBATCH --account pengel_spirit
#SBATCH --job-name readFiltering
#SBATCH --time 04:00:00
#SBATCH --mem 4G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/readFiltering.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/readFiltering.out

### Variables ###

username=gsartonl
workingDir=/users/${username}/spirit/${username}/assembly
yourfileID=10_S10
raw_data_directory=${workingDir}/00_rawReads

##!!!!!!!here you can adjust the parameters
##see how these bash variables are used in the script further bellow

MINIMUM_read_LENGTH=7000
min_mean_q=10
length_weight=10
target_bases=400000000

mkdir -p ${workingDir}/Oxford_Nanopore/01_filtlong_trimming
####--------------------------------------
##modules
####--------------------------------------

module load gcc/9.3.0
module load filtlong/0.2.0

####--------------------------------------
##Introduction to script
####--------------------------------------

start=$SECONDS
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${yourfileID}
echo "Step: Read filtering"

####--------------------------------------
##Filter reads
echo -e "1. First we filter the reads"



filtlong --min_length ${MINIMUM_read_LENGTH} \
    	 --min_mean_q ${min_mean_q} \
   		 --length_weight ${length_weight} \
    	 --target_bases ${target_bases}  \
    	${raw_data_directory}/${yourfileID}_filtered.fastq  > \
    	${workingDir}/Oxford_Nanopore/01_filtlong_trimming/${yourfileID}_filtered.1.fastq


### End ###
date +"END : %a %b %e %H:%M:%S %Z %Y"
duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
