#!/bin/bash
# USAGE: script.sh genomeindex-dir
#
# Using STAR for mapping NGS short-read data to genomes
# Here we do not do any mapping only removing the genome from mem

genomeindex=$1;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=./star_remove.${now}.stderr
outfile=./star_remove.${now}.stdout

STAR --outFileNamePrefix remove \
     --runThreadN 2 \
     --genomeLoad Remove \
     --genomeDir ${genomeindex} 2>> ${errfile} >> ${outfile};


