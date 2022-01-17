#!/bin/bash

echo "#####################################################"
echo " DRAM annotation : Starting time: $(date -u)"
echo "#####################################################"

# Setting for the new UTF-8 terminal support
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export PATH=/dcsrsoft/spack/external/dram/v1.2.4/bin:$PATH
module load gcc python
module load hmmer mmseqs2 prodigal infernal trnascan-se barrnap
# Define your variables
outDir=$1
checkMResult=$2;
Genomes=$3/contigs.fasta;
genomeID=$4;
echo  ${Genomes}
echo ${outDir}/${genomeID}
# Run DRAM
DRAM-setup.py prepare_databases --output_dir $3/../DRAM_data
DRAM.py annotate -i ${Genomes} -o ${outDir}/${genomeID} --min_contig_size 999 --threads 12 --verbose --checkm_quality ${checkMResult}
# cd ${outDir}/${genomeID}
# DRAM.py distill -i annotations.tsv -o distill --trna_path trnas.tsv --rrna_path rrnas.tsv
#

echo "#####################################################"
echo " Finishing time : $(date -u)"
echo "#####################################################"
