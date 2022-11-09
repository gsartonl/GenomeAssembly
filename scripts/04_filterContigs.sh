#!/bin/bash

### Filter the contigs based on length and K-mer coverage ##
contigs=$1;
coverage=$2;
length=$3;
outFile=$4

# store headers with $4= length and $6=coverage above specified values

echo "1) Getting contigs headers"
grep '>' ${contigs} | sed 's/_/ /g' | awk '{if($4>500 && $6>10){print $0}}' | sed 's/ /_/g' > $(dirname ${contigs})/HeaderFiltered.txt

echo "2) Extracting sequences"
# get the ids in an array - search pattern and if in ids copy until nex pattern
awk 'NR==FNR{ids[$0]; next} />/{f=($0 in ids)} f' $(dirname ${contigs})/HeaderFiltered.txt ${contigs} \
	> ${outFile};

echo "3) Removing temporary file"
rm $(dirname ${contigs})/HeaderFiltered.txt;

echo "### END ###"

### END ###
