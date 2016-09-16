#!/bin/bash
# USAGE: script.sh INDEX DIR-WITH-FASTQ SE/PE
# INDEX: SALMON-INDEX directory
# DIR-WITH-FASTQ: Will find all fastq in DIR-WITH-FASTQ including sub-directories
# SE/PE:  Either SE for single end or PE for paired-end.
# Looks for *.fastq.gz
# Will out put to directory ./quants

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
