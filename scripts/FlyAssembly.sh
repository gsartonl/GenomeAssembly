#!/bin/bash



#SBATCH --account pengel_spirit
#SBATCH --job-name Fly_assembly
#SBATCH --time 04:00:00
#SBATCH --mem 50G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/FlyAssembly.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/FlyAssembly.out

### Variables ###

username=gsartonl
workingDir=/users/${username}/spirit/${username}/assembly/Oxford_Nanopore
yourfileID=10_S10
GENOME_SIZE=2.1m

##!!!!!!!here you can adjust the parameters

rm -r ${workingDir}/02_assembly/${yourfileID}
mkdir -p ${workingDir}/02_assembly/${yourfileID}

start=$SECONDS


### modules ###


module load gcc/9.3.0
module load flye/2.8.3

date +"START : %a %b %e %H:%M:%S %Z %Y"
echo "Step: Genome assembly"


### Commands ###


flye --threads 8  --iterations 5 --genome-size ${GENOME_SIZE} \
      --nano-raw ${workingDir}/01_filtlong_trimming/${yourfileID}_filtered.1.fastq \
      --out-dir ${workingDir}/02_assembly/${yourfileID}


date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"


### END ###
