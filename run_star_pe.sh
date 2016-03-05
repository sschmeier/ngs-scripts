#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir

dir=$1;
outdir=$2;
genomeindex=$3;

now=$(date +"%Y-%d-%m_%H%M%S")

# set error and stdout file
errfile=${outdir}/star_pe.${now}.stderr
outfile=${outdir}/star_pe.${now}.stdout

outdir_2=${outdir}/star_pe.${now};
mkdir ${outdir_2};

# load genome into mem
STAR --outFileNamePrefix ${outdir_2} \
     --runThreadN 2 \
     --genomeLoad LoadAndExit \
     --genomeDir ${genomeindex} \
     2>> ${errfile} >> ${outfile};

for i in `ls ${dir}/*_1.fastq.trimmed.gz`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};

    nice STAR --readFilesIn ${i} $(echo ${i} | sed 's/_1/_2/g') \
         --outFileNamePrefix ${outdir_2}/$(basename ${i} | sed 's/_1.fastq.trimmed.gz//g') \
         --limitBAMsortRAM 20000000000 \
         --genomeLoad LoadAndKeep \
         --outFilterMultimapNmax 1 \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 2 \
         --alignIntronMax 1 \
         --alignEndsType EndToEnd \
         --genomeDir ${genomeindex} \
         --readFilesCommand zcat \
         2>> ${errfile} >> ${outfile};         

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
 
STAR --outFileNamePrefix ${outdir_2} \
     --runThreadN 2 \
     --genomeLoad Remove \
     --genomeDir ${genomeindex} \
     2>> ${errfile} >> ${outfile};

