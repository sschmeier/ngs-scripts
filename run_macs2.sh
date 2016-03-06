#!/bin/bash
# USAGE: script.sh in-dir out-dir genomesize

dir=$1;
outdir=$2;
genomesize=$3;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/macs2.${now}.stderr
outfile=${outdir}/macs2.${now}.stdout

outdir_2=${outdir}/macs2.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*.bam`; do

    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    macs2 callpeak -n ${outdir_2}/$(basename ${i} | sed 's/.bam//g') \
                   -t ${i} \
                   -f BAM \
                   -g ${genomesize} \
                   -q 0.05 \
                   2>> ${errfile} >> ${outfile}; 

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
 
gzip ${outdir_2}/*.narrowPeak
gzip ${outdir_2}/*.bed