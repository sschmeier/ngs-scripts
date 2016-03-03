#!/bin/bash

dir=$1;
outdir=$2;

now=$(date +"%Y-%d-%m_%H%M%S")

errfile=${outdir}/sickle_pe.${now}.stderr;
outfile=${outdir}/sickle_pe.${now}.stdout;

$outdir_2=${outdir_pre}/sickle_pe.${now}
mkdir ${outdir_2}

for i in `ls ${dir}/*_1*fastq*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};

    sickle pe \
-f ${i} \
-r $(echo ${i} | sed 's/_1/_2/g') \
-t sanger \
-g \
-o ${outdir_2}/$(basename ${i} | sed 's/.fastq.gz/.fastq.trimmed.gz/g') -p ${outdir_2}/$(basename ${i} | sed 's/_1/_2/g' | sed 's/.fastq.gz/.fastq.trimmed.gz/g') -s ${outdir_2}/$(basename ${i} | sed 's/_1/_singles/g' | sed 's/.fastq.gz/.fastq.trimmed.gz/g') \
2>> ${errfile} >> ${outfile};

    echo "-------------" >> ${errfile};
    echo "-------------" >> ${outfile};
done;
