#!/bin/bash
# USAGE: cat fastq | fastq-sort.sh > fastq.sorted

# input from stdin if no file is given 
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"
cat $input | paste -d',' - - - - | sort -k1,1 -S 8G | tr ',' '\n'
