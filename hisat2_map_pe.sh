#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-prefix removeMultiMap
# 
# hisat2 has no option to discard multimapping reads
# we report thus here all of them and need to filter afterwards with samtools
# output is a sorted bam-file.
#
# in-dir:             path to dir with fastq-files to map
#                     expects files of the form: *_1.fastq.gz, *_2.fastq.gz
# out-dir:            path to dir for results
# genomeindex-prefix: path+prefix to genomeindex
# removeMultiMap:     either "yes" or "no", if "yes" willl use samtools to remove multimappers
#

dir=$1;
outdir=$2;
genomeindex=$3;
rmmultimap=$4; 

now=$(date +"%Y-%m-%d_%H%M%S");

# set error and stdout file
errfile=${outdir}/hisat2_se.${now}.stderr;
outfile=${outdir}/hisat2_se.${now}.stdout;

outdir_2=${outdir}/hisat2_se.${now};
mkdir ${outdir_2};


for i in `ls ${dir}/*_1.fastq.gz`; do
    i2=$(echo ${i} | sed 's/_1/_2/g');
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    resultfile=${outdir_2}/$(basename $i | sed 's/.fastq.gz//g');
    
    # HISAT2 params:
    # --reorder: force SAM output order to match order of input reads
    # -q: query input files are FASTQ .fq/.fastq (default)
    # --all: report all alignments, very slow, MAPQ not meaningful
    # --no-spliced-alignment: disable spliced alignment
    
    if [ "$rmmultimap" = "yes" ]; then
        # use samtools to remove multimapper, sort, remove, duplicate reads
        # samtools view:
        #  -S = input format is auto-detected
        #  -u = uncompressed BAM output (implies -b)
        #  -h = include header in SAM output
        #  -q 1 = remove all reads with quality score 0, which are multimappers
        # samtools rmdup:
        #  -s = for single end reads
	nice hisat2 \
	    --threads 2 \
	    --reorder \
	    -q \
	    --phred33 \
	    --all \
	    -x ${genomeindex} \
	    --no-spliced-alignment \
	    -1 ${i} \
            -2 ${i2} 2>> ${errfile} | \
	    nice samtools view --threads 2 -Shu - 2>> ${errfile} | \
	    nice samtools sort --threads 2 - 2>> ${errfile} | \
	    nice samtools rmdup -s - ${resultfile}.sorted_nodup.bam >> ${outfile} 2>> ${errfile};
    else
	nice hisat2 \
 	    --threads 2 \
	    --reorder \
	    -q \
	    --phred33 \
	    --all \
	    -x ${genomeindex} \
	    --no-spliced-alignment \
	    -1 ${i} \
            -2 ${i2} 2>> ${errfile} | \
	    nice samtools view --threads 2 -bS - | \
	    nice samtools sort --threads 2 -o ${resultfile}.sorted.bam - 2>> ${errfile} >> ${outfile};
    fi
    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
