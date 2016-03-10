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

for i in `ls ${dir}/*_1.fastq.gz`; do
    fileprefix1=${outdir_2}/$(basename ${i})
    
    i2=$(echo ${i} | sed 's/_1/_2/g')
    fileprefix2=${outdir_2}/$(basename ${i2})

    echo "${i},${i2}" >> ${errfile};
    echo "${i},${i2}" >> ${outfile};

    # align
    nice bwa aln -t 4 ${genomeindex} ${i} -f ${fileprefix1}.sai >> ${outfile} 2>> ${errfile}
    nice bwa aln -t 4 ${genomeindex} ${i2} -f ${fileprefix2}.sai >> ${outfile} 2>> ${errfile}

    # remove multimapper, sort, remove,  duplicate reads
    nice bwa sampe ${genomeindex} ${fileprefix1}.sai ${fileprefix2}.sai ${i} ${i2} 2>> ${errfile} | \
	samtools view -Shu -q 1 - 2>> ${errfile} | \
	samtools sort - 2>> ${errfile} | \
	samtools rmdup - ${outdir_2}/$(basename ${i} | sed s/_1.fastq.gz/_sorted_nodup.bam/g) >> ${outfile} 2>> ${errfile}
    
    rm ${fileprefix1}.sai;
    rm ${fileprefix2}.sai;

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;


