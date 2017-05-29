#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir keepindex
#
# Using STAR for RNA-seq read alignment.
#
# in-dir:          path to directory with fastq.gz files to map
#                  expects files of the form: *_1.fastq.gz, *_2.fastq.gz
# out-dir:         path to directory for results 
# genomeindex-dir: path to directory with genomeindex
# keepindex:       yes/no; should the genomeindex be removed after run? 
#                  DEFAULT IS TO DROP INDEX = "no"
#

dir=$1;
outdir=$2;
genomeindex=$3;
keepindex=$4;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/star_pe.${now}.stderr
outfile=${outdir}/star_pe.${now}.stdout

outdir_2=${outdir}/star_pe.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*_1.fastq.gz`; do
    i2=$(echo ${i} | sed 's/_1/_2/g');
    echo ${i},${i2} >> ${errfile};
    echo ${i},${i2} >> ${outfile};
    
    # --limitBAMsortRAM NEEDS TO BE SET IF SHARED GENOME SHOULD BE USED
    # --outFilterMultimapNmax 1 = filter out multimappers
    # --genomeLoad LoadAndKeep = will load the genomindex specified 
    # if not already loaded and will keep it in MEM for subsequent runs
    nice STAR --readFilesIn ${i} ${i2} \
         --outFileNamePrefix ${outdir_2}/$(basename ${i} | sed 's/_1.fastq.gz/_/g') \
         --limitBAMsortRAM 20000000000 \
         --genomeLoad LoadAndKeep \
         --outFilterMultimapNmax 20 \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 4 \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outReadsUnmapped Fastx \
         --genomeDir ${genomeindex} \
         --readFilesCommand zcat 2>> ${errfile} >> ${outfile};      

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;

# Do we need to remove the genomeindex from MEM?
# Attention: If you use multiple runs at the same time, you might remove it too early.
if [ "$keepindex != "yes" ]; then 
    STAR --outFileNamePrefix ${outdir_2} \
         --runThreadN 2 \
         --genomeLoad Remove \
         --genomeDir ${genomeindex} 2>> ${errfile} >> ${outfile};
fi
