#!/bin/bash
# USAGE: script.sh genomeindex-dir
#
# Using STAR for mapping NGS short-read data to genomes
# Here we do not do any mapping only loading the genome into mem

genomeindex=$1;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=./star_load.${now}.stderr
outfile=./star_load.${now}.stdout

# load genome into mem
# MAKES ONLY SENSE IF --limitBAMsortRAM IS SET IN SUBSEQUENT STAR RUNS
STAR --outFileNamePrefix load \
     --runThreadN 2 \
     --genomeLoad LoadAndExit \
     --genomeDir ${genomeindex} 2>> ${errfile} >> ${outfile};


