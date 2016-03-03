#!/bin/bash

dir=$1;
outdir=$2;

now=$(date +"%Y-%d-%m_%H%M%S")

errfile=${outdir}/sickle_se.${now}.stderr;
outfile=${outdir}/sickle_se.${now}.stdout;

$outdir_2=${outdir_pre}/sickle_se.${now}
mkdir ${outdir_2}

for i in `ls ${dir}/*fastq*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    sickle se -f ${i} -t sanger -g -o ${outdir_2}/$(basename ${i} | sed 's/.fastq.gz/.fastq.trimmed.gz/g') 2>> ${errfile} >> ${outfile};
    echo "-------------" >> ${errfile};
    echo "-------------" >> ${outfile};
done;


