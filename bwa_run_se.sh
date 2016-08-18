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
    
    # align
    # -t 4 = use 4 threads
    # -f outfile
    nice bwa aln -t 4 ${genomeindex} ${i} -f ${fileprefix}.sai >> ${outfile} 2>> ${errfile}

    # remove multimapper, sort, remove,  duplicate reads
    # samtools view:
    #  -S = input format is auto-detected
    #  -u = uncompressed BAM output (implies -b)
    #  -h = include header in SAM output
    #  -q 1 = remove all reads with quality score 0, which are multimappers
    # samtools rmdup:
    #  -s = for single end reads
    nice bwa samse ${genomeindex} ${fileprefix}.sai ${i} 2>> ${errfile} | \
	samtools view -Shu -q 1 - 2>> ${errfile} | \
	samtools sort - 2>> ${errfile} | \
	samtools rmdup -s - ${outdir_2}/$(basename ${i} | sed s/.fastq.gz/_sorted_nodup.bam/g) >> ${outfile} 2>> ${errfile}
    
    rm ${fileprefix}.sai;

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;

