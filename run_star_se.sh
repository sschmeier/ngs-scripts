#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir

dir=$1;
outdir=$2;
genomeindex=$3;

now=$(date +"%Y-%d-%m_%H%M%S")

# set error and stdout file
errfile=${outdir}/star_se.${now}.stderr
outfile=${outdir}/star_se.${now}.stdout

outdir_2=${outdir}/star_se.${now};
mkdir ${outdir_2};

# load genome into mem
# MAKES ONLY SENSE IF --limitBAMsortRAM IS SET IN SUBSEQUENT RUNS
STAR --outFileNamePrefix ${outdir_2} \
     --runThreadN 2 \
     --genomeLoad LoadAndExit \
     --genomeDir ${genomeindex} \
     2>> ${errfile} >> ${outfile};

for i in `ls ${dir}/*.fastq.trimmed.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    # --limitBAMsortRAM NEEDS TO BE SET IF SHARED GENOME SHOULD BE USED
    nice STAR --outFileNamePrefix ${outdir_2}/$(basename ${i} | sed 's/.fastq.trimmed.gz//g') \
         --limitBAMsortRAM 20000000000 \
         --genomeLoad LoadAndKeep \
         --outFilterMultimapNmax 1 \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 2 \
         --alignIntronMax 1 \
         --alignEndsType EndToEnd \
         --genomeDir ${genomeindex} \
         --readFilesCommand zcat \
         --readFilesIn ${i} \
         2>> ${errfile} >> ${outfile}; 

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
 
STAR --outFileNamePrefix ${outdir_2} \
     --runThreadN 2 \
     --genomeLoad Remove \
     --genomeDir ${genomeindex} \
     2>> ${errfile} >> ${outfile};

