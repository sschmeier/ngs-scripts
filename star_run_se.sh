#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir yes/no
#
# Using STAR for mapping NGS short-read data to genomes, in a splice UNAWARE manner, similar to bwa, bowtie2, etc.
# Drop-in replacement for bwa, bowtie
# for RNAseq type of splice-aware mapping use differnt script.

dir=$1;
outdir=$2;
genomeindex=$3;
remove=$4

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/star_se.${now}.stderr
outfile=${outdir}/star_se.${now}.stdout

outdir_2=${outdir}/star_se.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*.fastq.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};
    
    # --limitBAMsortRAM NEEDS TO BE SET IF SHARED GENOME SHOULD BE USED
    # --outFilterMultimapNmax 1 = filter out multimappers
    # --alignIntronMax 1 = dont split up reads
    # --alignEndsType EndToEnd = dont split up reads
    # --genomeLoad LoadAndKeep = will load the genomindex specified if not already loaded and will keep it in MEM for subsequent runs
    nice STAR --outFileNamePrefix ${outdir_2}/$(basename ${i} | sed 's/.fastq.gz/_/g') \
         --limitBAMsortRAM 20000000000 \
         --genomeLoad LoadAndKeep \
         --outFilterMultimapNmax 1 \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 2 \
         --alignIntronMax 1 \
         --alignEndsType EndToEnd \
         --genomeDir ${genomeindex} \
         --readFilesCommand zcat \
         --readFilesIn ${i} 2>> ${errfile} >> ${outfile}; 

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;

# attention, if you use multiple runs at the same time, you might remove it too early. 
if [ "$remove" = "yes" ]; then
    STAR --outFileNamePrefix ${outdir_2} \
         --runThreadN 2 \
         --genomeLoad Remove \
         --genomeDir ${genomeindex} 2>> ${errfile} >> ${outfile};
fi

