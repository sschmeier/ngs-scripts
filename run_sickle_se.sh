#!/bin/bash

dir=$1;
outdir=$2;

# unique id
uuid=$(uuidgen);

errfile=${outdir}/sickle_se.${uuid}.stderr;
outfile=${outdir}/sickle_se.${uuid}.stdout;

$outdir_2=${outdir_pre}/sickle_se.${uuid}
mkdir ${outdir_2}

for i in `ls ${dir}/*fastq*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    sickle se -f ${i} -t sanger -g -o ${outdir_2}/$(basename ${i} | sed 's/.fastq.gz/.fastq.trimmed.gz/g') 2>> ${errfile} >> ${outfile};
    echo "-------------" >> ${errfile};
    echo "-------------" >> ${outfile};
done;


