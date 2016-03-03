#!/bin/bash
# USAGE: script.sh in-dir out-dir genomeindex-dir

dir=$1;
outdir=$2;
genomeindex=$3;

# unique id
uuid=$(uuidgen);

# set error and stdout file
errfile=${outdir}/star_se.${uuid}.stderr
outfile=${outdir}/star_se.${uuid}.stdout

outdir_2=${outdir}/star_se.${uuid};
mkdir ${outdir_2};

# load genome into mem
/mnt/DATA1/seb/bin/STAR_2.5.1b/STAR --runThreadN 2 --genomeLoad LoadAndExit --genomeDir ${genomeindex};

for i in `ls ${dir}/*.fastq.trimmed.gz*`; do
    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};

    /mnt/DATA1/seb/bin/STAR_2.5.1b/STAR --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --runThreadN 2 --alignIntronMax 1 --alignEndsType EndToEnd --genomeDir ${genomeindex} --readFilesCommand zcat --readFilesIn ${i} --outFileNamePrefix ${outdir_2}/$(basename ${i} | sed 's/.fastq.trimmed.gz//g' ) 2>> ${errfile} >> ${outfile}; 

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
 
/mnt/DATA1/seb/bin/STAR_2.5.1b/STAR --runThreadN 2 --genomeLoad Remove --genomeDir ${genomeindex};

