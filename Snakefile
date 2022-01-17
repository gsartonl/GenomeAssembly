#################################
# snakefile for genome assembly #
##          @gsartonl          ##
##          06/01/2022         ##
#################################

## Import config file ##
configfile: "config.yaml"

### Loading modules ###
import os
from datetime import datetime
from pathlib import Path


##
LIBRARY=config['SAMPLES']
#print(LIBRARY)
## SetUp paths ##

workDir = config['directories']['workDir']
pathRaw = workDir + config['directories']['rawReadsDir']
pathFastQC = workDir + config['directories']['fastQC']
pathTrimmed = workDir + config['directories']['trimmedReads']
pathAssembly = workDir + config['directories']['assembly']
pathAnnotation = workDir + config['directories']['annotation']
pathCheckM = workDir + config['directories']['checkM']
pathCoverage = workDir + config['directories']['readCoverage']
coverage=config['FilterContigs']['coverage']
length=config['FilterContigs']['length']

#print(pathTrimmed)

### Functions ###
def GetPaired_input(wildcards) :
   fwd=pathRaw+"/"+wildcards.LIBRARY+"_L"+wildcards.LANE+"_R1_001_"+config['metadata'][wildcards.LIBRARY][int(wildcards.LANE)]['Hash_FWD']+".fastq.gz"
   rev=pathRaw+"/"+wildcards.LIBRARY+"_L"+wildcards.LANE+"_R2_001_"+config['metadata'][wildcards.LIBRARY][int(wildcards.LANE)]['Hash_REV']+".fastq.gz"
  # print([fwd , rev])
   return([fwd , rev])

# Novogene formated
def GetPaired_input(wildcards) :
    #print(wildcards.LIBRARY)
    fwd=pathRaw+"/"+wildcards.LIBRARY+ "_" +config['LIBRARY'][wildcards.LIBRARY]['Hash_FWD']+ "_L"+config['LIBRARY'][wildcards.LIBRARY]['Lane']+"_1" + config['suffix']
    rev=pathRaw+"/"+wildcards.LIBRARY+ "_" +config['LIBRARY'][wildcards.LIBRARY]['Hash_REV']+ "_L"+config['LIBRARY'][wildcards.LIBRARY]['Lane']+"_2" + config['suffix']
    #print([fwd , rev])
    return([fwd , rev])

def GetGenomeID(path):
    fileName=path.strsplit('/')[1]
    ID=fileName.strsplit('_')[0]
    return(ID)

## rules ##
#print(LIBRARY)
rule all :
    input :
        pairedFWD=expand(pathTrimmed + "/{LIBRARY}_R1_paired.fastq.gz", LIBRARY="ESL0744"),
        outDir=expand(pathFastQC + "/{LIBRARY}", LIBRARY=LIBRARY),
        outDirAssembly=expand(pathAssembly + "/{LIBRARY}", LIBRARY=LIBRARY),
        #checkM=expand(pathCheckM + "/{LIBRARY}/", LIBRARY=LIBRARY),
        tsv=expand(pathCheckM + "/{LIBRARY}/CheckM_QC_f2_{LIBRARY}.tsv", LIBRARY=LIBRARY),
        outDirDRAM= expand(pathAnnotation + "/{LIBRARY}", LIBRARY=LIBRARY),
        depth=expand(pathCoverage + "/{LIBRARY}_depth.txt", LIBRARY=LIBRARY)

rule fastQC :
    input:
        GetPaired_input
    output:
        outDir=directory(pathFastQC + "/{LIBRARY}")

    resources:
        mem_mb=20000,
        account='pengel_spirit'
    threads: 4
    log:
        "logs/01_fastcq/fastqc_{LIBRARY}"
    message : "Runing FastQC"
    conda:
        "envs/01_fastqc.yaml"
    shell:
        "scripts/fastqc.sh {input} {output.outDir}"


rule readTrimming:
    input:
        GetPaired_input
    output:
        pairedFWD=pathTrimmed + "/{LIBRARY}_R1_paired.fastq.gz",
        pairedREV=pathTrimmed + "/{LIBRARY}_R2_paired.fastq.gz",
        unpairedFWD=pathTrimmed + "/unpaired/{LIBRARY}_R1_unpaired.fastq.gz",
        unpairedREV=pathTrimmed + "/unpaired/{LIBRARY}_R2_unpaired.fastq.gz"

    params:
        outDir = pathTrimmed,
        genomeID = '{LIBRARY}',
        adaptors = ''.join(config['illuminaAdaptors'])
    resources:
        runtime_s="3600", # memory in secs
        mem_mb=20000, # double memory if job fails
        account='pengel_spirit'
    threads: 8
    message:"Runing trimmomatic"
    log:
        "logs/02_trimmed/ReadTrimming_{LIBRARY}"
    #     "logs/02_trimmed/ReadTrimming_{LIBRARY}.cluster.e"
    conda:
        "envs/02_trimmomatic.yaml"
    shell:
        "scripts/02_trimmomatic.sh {input} {params.outDir} {params.genomeID} {params.adaptors}"

rule assembly:
    input:
        FWD=rules.readTrimming.output.pairedFWD,
        REV=rules.readTrimming.output.pairedREV
    output:
        outDir=directory(pathAssembly + "/{LIBRARY}"),
        contigs=pathAssembly + "/{LIBRARY}/contigs.fasta"
    params:
        inDir=pathTrimmed + "/merged",
        genomeID="{LIBRARY}"
    resources:
        account="pengel_spirit",
        runtime_s="86400",
        mem_mb=700000
        #mem_mb=lambda wildcards, attempt:700000+attempt*10000
    threads: 16
    message: "SPADes assembly"
    log:
        "logs/03_assembly/assembly_{LIBRARY}"
    conda:
        "envs/03_SPAdes.yaml"
    shell:
        "scripts/03_spades.sh {input.FWD} {input.REV} {output.outDir} {threads}"

rule FilterContigs :
    input:
        assembly=rules.assembly.output.contigs
    output:
        filteredContigs= pathAssembly + "/{LIBRARY}/contigs_filtered.fasta"
    params:
       coverage=coverage,
       length=length
    resources:
        account='pengel_spirit',
        runtime_s="600",
        mem_mb=1500
    log:
        "logs/04_FilterContigs/FilterContigs_{LIBRARY}"
    shell:
        "scripts/04_filterContigs.sh {input.assembly} {params.coverage} {params.length} {output.filteredContigs}"


rule checkM:
    input:
        assembly=rules.FilterContigs.output.filteredContigs
    output:
        #outDir=directory(pathCheckM + "/{LIBRARY}/")
        tsv=pathCheckM + "/{LIBRARY}/CheckM_QC_{LIBRARY}.tsv",
        tsvf2=pathCheckM + "/{LIBRARY}/CheckM_QC_f2_{LIBRARY}.tsv"

    params:
        workDir=workDir,
        outDir=pathCheckM + "/{LIBRARY}",
        genomeID="{LIBRARY}"
    resources:
        account='pengel_spirit',
        runtime_s="86400",
        mem_mb=150000
    log:
        "logs/05_checkM/checkM_{LIBRARY}"
    conda:
        "envs/05_checkM.yaml"
    shell:
        "scripts/05_checkM.sh {params.workDir} {input.assembly} {params.outDir} {params.genomeID}"

# works only on the cluster
rule dram:
    input:
        assembly=rules.assembly.output.outDir,
        tsvf2=rules.checkM.output.tsv
    output:
       outDir= directory(pathAnnotation + "/{LIBRARY}")
    params:
      outDir=pathAnnotation,
      genomeID="{LIBRARY}"
    resources:
        account='pengel_spirit',
        runtime_s="1800",
        mem_mb=150000
    threads: 16
    log : "logs/06_DRAM/DRAM_{LIBRARY}"
    message: "DRAM annotation"
    #conda:
    #    "envs/05_dram.yaml"
    shell:
        "scripts/06_dram.sh {params.outDir} {input.tsvf2} {input.assembly} {params.genomeID}"

rule readMapping:
    input:
        FWD=rules.readTrimming.output.pairedFWD,
        REV=rules.readTrimming.output.pairedREV,
        assembly=rules.FilterContigs.output.filteredContigs
    output:
      sam=pathCoverage + "/{LIBRARY}_mapping.sam",
      bam=pathCoverage + "/{LIBRARY}_mapping.bam",
      sorted=pathCoverage + "/{LIBRARY}_mapping_sorted.bam",
      depth=pathCoverage + "/{LIBRARY}_depth.txt"

    params:
      genomeID="{LIBRARY}"
    resources:
        account='pengel_spirit',
        runtime_s="1800",
        mem_mb=150000 #150G
    threads: 16
    log : "logs/QC/QC_{LIBRARY}"
    conda:
      "envs/07_readCoverage.yaml"
    shell:
      "scripts/07_readMapping.sh {input.assembly} {input.FWD} {input.REV} {params.genomeID} {output.sam} {output.bam} {output.sorted} {output.depth}"

# rule plotCoverage:
#     input:
#         depth=rules.readCoverage.ouput.depth
#     output:
#         depthPlot=pathCoverage + "/{LIBRARY}_coveragePlot.png"
#     params :
#       coverage=config['FitterContigs']['coverage']
#       length=config['FitterContigs']['length']
#     resources:
#         account='pengel_spirit',
#         runtime_s="600",
#         mem_mb=1500 #1.5G
#     log : "logs/QC/QC_{LIBRARY}"
#     conda:
#       "envs/08_plotCoverage.yaml"
#     shell:
#       "scripts/08_plotCoverage.sh {input.assembly} {params.coverage} {params.length} {output.filteredContigs}"
