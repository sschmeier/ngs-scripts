#!/bin/bash
dir=$1;
outdir=$2;

now=$(date +"%Y-%d-%m_%H%M%S")

errfile=${outdir}/fastqc.${now}.stderr;
outfile=${outdir}/fastqc.${now}.stdout;

outdir_2=${outdir}/fastqc.${now};
mkdir ${outdir_2};

fastqc -q -o ${outdir_2} ${dir}/*fastq* 2>> ${errfile} >> ${outfile};
rm ${outdir_2}/*.html;





