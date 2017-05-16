# Notes on how to perform a salmon-based tx analysis

## Transcript gene mapping
Get ensembl gene to trasncript mapping in accordaance to the same genome version used for salmon. For example:

```bash
# download ensemble
nohup lftp -c "open ftp://ftp.ensembl.org/pub/release-85/gtf/homo_sapiens; mirror ./ ./Homo_sapiens-release-85/"

# extract
zcat ./Homo_sapiens-release-85/Homo_sapiens.GRCh38.85.gtf.gz | cut -f 9 | seb.FILE.grep - '(ENSG.+?)".+(ENST.+?)".+transcript_version\s"(.+?)".+gene_name\s"(.+?)"' | sort | uniq | awk -F '\t' '{print $2"."$3"\t"$1"\t"$4}' | sort -k3,3 > tx_gene_map.txt
```

## Salmon indexing
One might want to consider all cdna + all ncrna. These fasta-files need to be combined in such a case. However, here an example for cdna.

```bash
# get fasta
nohup lftp -c "open ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cdna; mirror ./ ./Homo_sapiens-release-85-fasta/"
gzip -d ./Homo_sapiens-release-85-fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz

# create salmon index
nice salmon index -t ./Homo_sapiens-release-85-fasta/Homo_sapiens.GRCh38.cdna.all.fa -i Homo_sapiens.GRCh38.cdna.all.idx -p 4
```

## Salmon paired-end pseudo-mapping
Use the file `salmon_quant_multiple.sh` to perform the mapping for all files in directory `mates/`. This directory contains all samples that need to be mapped. Two files per sample (ending in `_1.fastq.gz` and `_2.fastq.gz`), i.e. belonging to a paired-end reads. 

```bash
bash salmon_quant_multiple.sh Homo_sapiens.GRCh38.cdna.all.idx mates/ PE tx_gene_map.txt
```

## TPM quantification of genes / transcripts
Two options are recommended. 1/ Use tximport `lengthscaled` method, 2/ use edgeR method. Log-transformation of data is optional.

```bash
# Do not do any log transform here.
Rscript salmon_get_genes_TPM.R lengthScaledTPM no | gzip > genes_salmon_tximport_TPM-lengthscaled.txt.gz
Rscript salmon_edgeR_genes_TPM.R no | gzip > genes_salmon_edgeR_TPM.txt.gz

# for transcipts
Rscript salmon_get_tx_TPM.R lengthScaledTPM no | gzip > tx_salmon_tximport_TPM-lengthscaled.txt.gz
Rscript salmon_edgeR_tx_TPM.R no | gzip > tx_salmon_edgeR_TPM.txt.gz
```
