#!/bin/bash -l

#Prerequisites: install methylseq pipeline in current directory 
export NXF_HOME="Path to methylseq pipeline directory"

NXF_OPTS='-Xms1g -Xmx4g'


gtf="../annotation/vcard.gtf"
fasta="../reference/GCA_905220365.1_ilVanCard2.1_genomic_chroms.fasta"

nextflow run nf-core/rnaseq -name "rnaseq" --input sample_sheet.csv -profile uppmax --project snic2022-5-34 --max_cpus 20 --max_memory 128GB --fasta $fasta --gtf $gtf --outdir ./results2 --skip_preseq --skip_stringtie --clip_r1 9 --clip_r2 9 --three_prime_clip_r1 7 --three_prime_clip_r2 7
