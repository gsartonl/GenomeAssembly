########## Welcome to the lab genome assembler ##########
###					 	  ###
### 		@gsartonl 			  ###
###                                               ###



There are a few setup steps before runing the pipeline.

## 1) Prepare the configuration file ##

First, open the config.yaml in the `setup` directory
-> change the path for the main output directory under the variable `workDir`.

Second, adapt the raw file `suffix` if needeed.

## 2) Add raw files ##

Create a directory with the SAME name as the one you put in the config file. 

-> in the `<outDir>`, create a sub-directory with the name `00_rawReads`
-> Add the raw sequencing files in the sub directory.

## 3) Activate the snakemake environment ##

Open a new terminal window
-> `cd ~/GenomeAssembly/`

Activate the snakemake environment
-> `conda activate snakemake` 

## 4) Prepare the database ##

Paste the following command in the terminal to set-up the database replace `<outDir>` by
the name you selected previously
-> `python scripts/00_SetUpDB.py ~/GenomeAssembly/ ./<outDir>/00_rawReads/`

## 5) Run snakemake pipeline ##

Copy and paste the following command line in the terminal
-> `snakemake --use-conda --core 4 all`

The pipeline takes a bit of time as the computer is not the most powerfull. The more genomes
you have the longer it will take (it will not take days but 5 to 10 mins per genome - if possible
some process are executed in parallel)

## 6) Get your results ##

To find your results, open Finder and click on `GenomeAssembly` in the favourites.

Here you will find : 
 - 00_rawReads : raw reads 
 - 01_FastQC : FastQC reports for each sample
 - 02_trimmed : 
 - 03_assembly : in a directory per sample : SPAdes assembly results, including contigs.fasta and contigs_filtered.fasta files
		Contigs with a coverage <10 and/or length < 500 are removed in the contigs_filtered.fasta file
 - 04_annotation :  in a directory per sample : Prokka annotation results. File names and locus tags are the sample name



If you have any troubles, send me an email : garance.sarton-loheac@unil.ch or contact me via slack
