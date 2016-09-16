#!/bin/bash
# USAGE: scipt.sh DIR 

DIR="$1"
for SRA in `find $DIR -name '*.sra'`; do
    echo "Processing:" $SRA
    DIR=$(dirname "${SRA}")
    fastq-dump --outdir ${DIR} --split-files ${SRA}
    gzip ${DIR}/*.fastq
done
