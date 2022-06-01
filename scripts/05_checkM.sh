#!/bin/bash --login


echo "#####################################################"
echo $id "Starting time: $(date -u)"
echo "#####################################################"


################# Setting up CheckM #################
checkm data setRoot "/Users/sbrochet/GenomeAssembly/setup/checkm_databases"

checkm test ./checkm_test_results
# ### Create directories and define variables
workdir=$1
assembly=$2
outDir=$3;
genomeID=$4;

plt="${outDir}/plots"

mkdir ${outDir}
mkdir ${plt}

echo -e "\n\nYour annotation working directory is ${workdir}\n"
echo -e "Your CheckM directory should be here: ${out} \n"
echo -e "Your CheckM plots will be saved here: ${plt}\n \n"

################# Running CheckM #################
# Run CheckM with the lineage workflow (option --genes means this are AA instead of NT)
checkm lineage_wf -t 24 --tab_table -f ${outDir}/CheckM_QC_${LIBRARY}.tsv -x fasta --aai_strain 0.95 ${assembly} ${outDir}
# Generate the extended table for DRAM
# checkm qa -t 24 -o 2 --tab_table -f ${outDir}/CheckM_QC_f2_${LIBRARY}.tsv ${outDir}/lineage.ms ${outDir}
# Make Plots!!!
#checkm gc_plot ${workdir} ${plt} -x fasta 95
#checkm coding_plot -x fasta ${out} ${bin} ${plt} 95 ${out}
#    checkm nx_plot ${bin} ${plt} -x fasta ${out}
#checkm marker_plot ${out} ${bin} ${plt}

#rm -rf /scratch/jgianott/sage/SAGE2021_22/${user}/DRAMkefir/checkm_databases/

echo "#####################################################"
echo $id " Finishing time : $(date -u)"
echo "#####################################################"
