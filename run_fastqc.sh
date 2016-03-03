#!/bin/bash
dir=$1;
outdir=$2;

# unique id
uuid=$(uuidgen);

errfile=${outdir}/fastqc.${uuid}.stderr;
outfile=${outdir}/fastqc.${uuid}.stdout;

outdir_2=${outdir}/fastqc.${uuid};
mkdir ${outdir_2};

fastqc -q -o ${outdir_2} ${dir}/*fastq* 2>> ${errfile} >> ${outfile};
rm ${outdir_2}/*.html;





