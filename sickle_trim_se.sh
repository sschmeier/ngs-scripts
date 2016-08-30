#!/bin/bash

dir=$1;
outdir=$2;

now=$(date +"%Y-%m-%d_%H%M%S")

errfile=${outdir}/sickle_se.${now}.stderr;
outfile=${outdir}/sickle_se.${now}.stdout;

outdir2=${outdir}/sickle_se.${now};
mkdir ${outdir2};

for i in `ls ${dir}/*fastq*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    # sickle:
    #  -f, --fastq-file, Input fastq file (required)
    #  -t, --qual-type,  = sanger
    #  -o, --output-file, Output trimmed fastq file (required)
    #  -g, --gzip-output, Output gzipped files.
    sickle se -f ${i} -t sanger -g -o "${outdir2}/$(basename ${i} | sed 's/.fastq.gz/.trimmed.fastq.gz/g')" 2>> ${errfile} >> ${outfile};
    echo "-------------" >> ${errfile};
    echo "-------------" >> ${outfile};
done;


