#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir
# 
# hisat2 has no option to discard multimapping reads
# we report thus here all of them and need to filter afterwards
# with samtools


dir=$1;
outdir=$2;
genomeindex=$3;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/hisat2_se.${now}.stderr
outfile=${outdir}/hisat2_se.${now}.stdout

outdir_2=${outdir}/hisat2_se.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*.fastq.trimmed.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    resultfile=${outdir_2}/$(basename $i | sed 's/.fastq.trimmed.gz/.sam/g')
    
    nice hisat2 \
     --threads 2 \
     --reorder \
     -q \
     --phred33 \
     --end-to-end \
     --all \
     -x ${genomeindex} \
     -U ${i} \
     -S ${resultfile} \
     2>> ${errfile} >> ${outfile}; 
    
    #nice samtools view -bS ${resultfile} > ${resultfile}.bam
    #rm ${resultfile}
    
    #samtools sort ${resultfile}.bam ${resultfile}.sorted.bam
    #rm ${resultfile}.bam

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
