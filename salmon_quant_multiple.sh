#!/bin/bash
# USAGE: script.sh INDEX DIR-WITH-FASTQ SE/PE
#
# INDEX: 
#   SALMON-INDEX directory
#   Create index withm e.g.
#   salmon index -t ../../TX/Mus_musculus-release-85/cdna.ncrna/Mus_musculus.GRCm38.cdna.ncrna.fa -i Mus_musculus.GRCm38.cdna.ncrna.idx -p 4
#
# DIR-WITH-FASTQ: Will find all fastq-files in DIR-WITH-FASTQ including sub-directories
#
# SE/PE: Either SE for single-end/unpaired or PE for paired-end reads (expects _1.fastq.gz and _2.fastq.gz files) .
#
# Looks for *.fastq.gz
# Will write output to directory ./quants
#

INDEX=$1
DIR=$2
TYPE=$3
# create base out dir
mkdir -p ./quants

for FQ in $(find $DIR -name '*.fastq.gz'); do
    if [[ $TYPE == "PE" ]]; then	
        # skip if _2.fastq in filename as these are processed with _1.fastq files
        if [[ $FQ == *"_2.fastq"* ]]; then
	    echo "Skipping as processed with read 1:" `basename $FQ`;
	    continue
        fi

        FQ1=${FQ}
        FQ2=${FQ1/_1/_2}

	TEMP1=$(basename ${FQ1}) 
	TEMP=${TEMP1%%.*} # remove all after . 
	DIR=./quants/${TEMP/_1/}	
	mkdir -p ${DIR}

	echo PE processing: ${FQ1}, ${FQ2}, Results in: $DIR
	# RUN SALMON
	# --useVBOpt: VB optimizer rather than the standard EM optimizer. 
        # tends to give more accurate results than the standard EM algorithm.
	salmon quant -l A \
	             -i ${INDEX} \
                     -1 <(gunzip -c ${FQ1}) \
                     -2 <(gunzip -c ${FQ2}) \
                     -o ${DIR} \
                     -p 4 \
	             --useVBOpt \
                     --numBootstraps 30 > ${DIR}/salmon.stdout 2> ${DIR}/salmon.stderr        
    
    else # SE
	TEMP1=$(basename ${FQ}) 
	TEMP=${TEMP1%%.*} # remove all after . 
	DIR=./quants/${TEMP/_1/}	
	mkdir -p ${DIR}

	echo SE processing: ${FQ}, Results in: $DIR
	# RUN SALMON
	# --useVBOpt: VB optimizer rather than the standard EM optimizer. 
        # tends to give more accurate results than the standard EM algorithm.
	salmon quant -l A \
	             -i ${INDEX} \
                     -r <(gunzip -c ${FQ}) \
                     -o ${DIR} \
                     -p 4 \
                     --useVBOpt \
                     --numBootstraps 30 > ${DIR}/salmon.stdout 2> ${DIR}/salmon.stderr        
    fi
done
