#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-prefix rm-multi-dup
#
# Using bwa for mapping NGS short-read data to genomes, 
# 
# in-dir:             path to directory with fastq.gz files to map
# out-dir:            path to directory for results 
# genomeindex-prefix: path to genomeindex-prefix and index
# rm-multi-dup:       either "yes" or "no", if "yes" will use samtools to remove multimapper
#                     and duplicates
#

indir=$1;
outdir=$2;
genomeindex=$3;
rmmultimap=$4; 

now=$(date +"%Y-%m-%d_%H%M%S");

# set error and stdout file
errfile=${outdir}/bwa_se.${now}.stderr;
outfile=${outdir}/bwa_se.${now}.stdout;

outdir_2=${outdir}/bwa_se.${now};
mkdir ${outdir_2};

for i in `ls ${indir}/*.fastq.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    fileprefix=${outdir_2}/$(basename $i | sed 's/.fastq.gz//g');
    
    if [ "$rmmultimap" = "yes" ]; then
        # remove multimapper, sort, remove, duplicate reads
        # samtools view:
        #  -S   : input format is auto-detected
        #  -u   : uncompressed BAM output (implies -b)
        #  -h   : include header in SAM output
        #  -q 1 : remove all reads with quality score 0, which are multimappers
        #  -b   : output format bam 
        # samtools rmdup:
        #  -s   : for single end reads
       
        nice bwa mem -t 4 ${genomeindex} ${i} 2>> ${errfile} | \
	        nice samtools view -Shu -q 1 - 2>> ${errfile} | \
	        nice samtools sort --threads 2 - 2>> ${errfile} | \
	        nice samtools rmdup -s - ${fileprefix}.sorted_nodup.bam >> ${outfile} 2>> ${errfile};
        
    else
        nice bwa mem -t 4 ${genomeindex} ${i} 2>> ${errfile} | \
            nice samtools view --threads 2 -bS - 2>>  ${errfile} | \
	        nice samtools sort --threads 2 -o ${resultfile}.sorted.bam - 2>> ${errfile} >> ${outfile};
    fi
    
    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;


