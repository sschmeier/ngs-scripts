#!/bin/bash
# USAGE: script.sh in-dir out-dir genome.fa

dir=$1;
outdir=$2;
genomeindex=$3;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/bwa_se.${now}.stderr
outfile=${outdir}/bwa_se.${now}.stdout

outdir_2=${outdir}/bwa_se.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*.fastq.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    fileprefix=${outdir_2}/$(basename ${i})

    nice bwa aln -t 4 ${genomeindex} ${i} > ${fileprefix}.sai 
    nice bwa samse ${genomeindex} ${fileprefix}.sai ${i} | \
	samtools view -Shu -q 1 - | \
	samtools sort - - | \
	samtools rmdup -s - - | \
	> ${fileprefix}_sorted_nodup.bam 

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;


