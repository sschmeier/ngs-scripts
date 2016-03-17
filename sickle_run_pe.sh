#!/bin/bash

dir=$1;
outdir=$2;

now=$(date +"%Y-%m-%d_%H%M%S")

errfile=${outdir}/sickle_pe.${now}.stderr;
outfile=${outdir}/sickle_pe.${now}.stdout;

outdir_2=${outdir}/sickle_pe.${now}
mkdir ${outdir_2}

for i in `ls ${dir}/*_1*fastq.gz`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};

    # sickle:
    #  -f, --pe-file1, Input paired-end forward fastq file (Input files must have same number of records)
    #  -r, --pe-file2, Input paired-end reverse fastq file
    #  -o, --output-pe1, Output trimmed forward fastq file
    #  -p, --output-pe2, Output trimmed reverse fastq file. Must use -s option.
    #  -t, --qual-type,  = sanger
    #  -g, --gzip-output, Output gzipped files.

    sickle pe \
-f ${i} \
-r $(echo ${i} | sed 's/_1/_2/g') \
-t sanger \
-g \
-o ${outdir_2}/$(basename ${i} | sed 's/.fastq.gz/.trimmed.fastq.gz/g') -p ${outdir_2}/$(basename ${i} | sed 's/_1/_2/g' | sed 's/.fastq.gz/.trimmed.fastq.gz/g') -s ${outdir_2}/$(basename ${i} | sed 's/_1/_singles/g' | sed 's/.fastq.gz/.trimmed.fastq.gz/g') 2>> ${errfile} >> ${outfile};

    echo "-------------" >> ${errfile};
    echo "-------------" >> ${outfile};
done;
