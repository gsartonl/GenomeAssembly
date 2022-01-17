#!/bin/bash

#SBATCH --account pengel_spirit
#SBATCH --job-name racon_polishing
#SBATCH --time 04:00:00
#SBATCH --mem 8G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/raconPolishing.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/raconPolishing.out


### Variables ###

username=gsartonl
yourfileID=10_S10

workingDir=/users/${username}/spirit/${username}/assembly/Oxford_Nanopore
outDir=${workingDir}/03_Polishing_ONT/${yourfileID}
reference_fasta=/users/${username}/spirit/${username}/assembly/test.fa
#reference_fasta=${workingDir}/02_assembly/10_S10/assembly.fasta
ont_reads=${workingDir}/01_filtlong_trimming/${yourfileID}_filtered.1.fastq

####--------------------------------------
### modules ###


module load gcc/9.3.0
module load graphmap/0.6.4
module load racon/2021-06-02

### Commands ###


start=$SECONDS
echo "Step: Racon polishing"
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${yourfileID}


### Priming the variables used and starting the for-loop ###
 
overall_polishing_counter=1
rm -r ${outDir}

for ONT_polishing_counter in $(echo "first second") ;
  do 
  echo -e "----"${ONT_polishing_counter}" round of Racon polishing"

### graphmap mapping ###

##only run this line if you do not have a fastq file yet. 
#gunzip -c ${personal_home_directory}/01_data/After_filtlong_trimming/${sample_name}_filtered.fastq.gz > ${personal_home_directory}/01_data/After_filtlong_trimming/${sample_name}_filtered.fastq 

#mkdir -p ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/{Graphmap_ONT,assembly}

#graphmap align --rebuild-index --circular  \
#  -r $reference_fasta \
#  -d $ont_reads \
#  -o ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam 

###===========================
##racon polishing
echo -e "-------2. Second polish with Racon"
echo ${ont_reads}
echo ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam
echo ${reference_fasta}
###===========================   

racon ${ont_reads} \
  ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}_Racon_ONT2FinalAssembly_sorted.sam \
  ${reference_fasta} > ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${yourfileID}_${ONT_polishing_counter}_Racon_Assembly.fasta

###===========================
##reset the parameters
echo -e "-------3. Resetting the parameters"
###===========================   

reference_fasta=$(echo ${outDir}/0${overall_polishing_counter}_${ONT_polishing_counter}_Racon_graphmap/assembly/${yourfileID}_${ONT_polishing_counter}_Racon_Assembly.fasta)
overall_polishing_counter=$((overall_polishing_counter+1)) 

done #polishing round


####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
