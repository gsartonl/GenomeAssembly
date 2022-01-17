#!/bin/bash

#SBATCH --account pengel_spirit

#SBATCH --job-name pilon_polishing
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --mem 12G
#SBATCH --time 03:30:00
#SBATCH --error /users/gsartonl/spirit/gsartonl/assembly/logs/pilonPolishing.err
#SBATCH --output /users/gsartonl/spirit/gsartonl/assembly/logs/pilonPolishing.out
####--------------------------------------${sample}
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------

username=gsartonl
yourfileID=10_S10
threads=6
workingDir=/users/${username}/spirit/${username}/assembly
Overall_output_directory=${workingDir}/Oxford_Nanopore/04_Polishing_Illumina/${yourfileID}
reference_fasta=${workingDir}/Oxford_Nanopore/03_Polishing_ONT/${yourfileID}/02_second_Racon_graphmap/assembly/10_S10_second_Racon_Assembly.fasta
#/users/gsartonl/spirit/gsartonl/assembly/Oxford_Nanopore/03_Polishing_ONT/10_S10/02_second_Racon_graphmap/assembly/_second_Racon_Assembly.fasta

raw_data_directory_illumina=${workingDir}/00_rawReads

####--------------------------------------
##modules
####--------------------------------------

module load gcc/9.3.0
module load pilon/1.24
module load bowtie2/2.4.2
module load samtools/1.12

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: Pilon polishing"

#done
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${yourfileID}

###===========================
##Priming the variables used and starting the for-loop
echo -e "-------0. Priming the variables used and starting the for-loop"
###===========================   

overall_polishing_counter=3

rm -r ${Overall_output_directory}

for Illumina_polishing_counter in $(echo "first second third")
do 

echo -e "----"${Illumina_polishing_counter}" round of pilon polishing"

###===========================
##bowtie2 mapping
echo -e "-------1. First we map with bowtie2"
###===========================      
mkdir -p $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/{bowtie2,assembly}

bowtie2-build --quiet $reference_fasta $reference_fasta

bowtie2 -1 ${raw_data_directory_illumina}/*_${yourfileID}_*_R1*.fastq.gz \
        -2 ${raw_data_directory_illumina}/*_${yourfileID}_*_R2*.fastq.gz \
        -x  $reference_fasta \
        --threads ${threads} \
        --local --very-sensitive-local | samtools sort -O BAM -o $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam - 
     
###--------------------
##index bam file
###--------------------     
samtools index $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam

###===========================
##pilon polishing
echo -e "-------2. Second polish with pilon"
###===========================   

pilon --threads ${threads} --genome $reference_fasta \
  --frags $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam \
  --output ${yourfileID}_${Illumina_polishing_counter}_Pilon_Assembly \
  --outdir $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/assembly/ \
  --changes

###===========================
##reset the parameters
echo -e "-------3. Resetting the parameters"
###===========================   

reference_fasta=$(echo $Overall_output_directory/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/assembly/${yourfileID}_${Illumina_polishing_counter}_Pilon_Assembly.fasta)
echo $reference_fasta
overall_polishing_counter=$((overall_polishing_counter+1)) 

done #polishing round

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
