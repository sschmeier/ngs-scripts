#!/bin/bash
dir=$1;
outdir=$2;

# unique id
uuid=$(uuidgen);

errfile=${outdir}/fastqc.${uuid}.stderr;
outfile=${outdir}/fastqc.${uuid}.stdout;

outdir_2=${outdir}/fastqc.${uuid};
mkdir ${outdir_2};

fastqc -q --extract -o ${outdir_2} ${dir}/*fastq* 2>> ${errfile} >> ${outfile};
rm ${outdir_2}/*.html;
rm ${outdir_2}/*.zip;

# look for warnings and fails
numlines_fail=$(find ${outdir_2} -name summary.txt | xargs -I pattern egrep '^FAIL.*Overrepresented' pattern | wc -l);
numlines_warn=$(find ${outdir_2} -name summary.txt | xargs -I pattern egrep '^WARN.*Overrepresented' pattern | wc -l);
echo "FAIL OVERREPRESENTED: ${numlines_fail}" >> ${outfile};
echo "WARN OVERREPRESENTED: ${numlines_fail}" >> ${outfile};

numlines_fail=$(find ${outdir_2} -name summary.txt | xargs -I pattern egrep '^FAIL.*Adapter' pattern | wc -l);
numlines_warn=$(find ${outdir_2} -name summary.txt | xargs -I pattern egrep '^WARN.*Adapter' pattern | wc -l);
echo "FAIL ADAPTER: ${numlines_fail}" >> ${outfile};
echo "WARN ADAPTER: ${numlines_fail}" >> ${outfile};

echo "DONE" >> ${outfile};



