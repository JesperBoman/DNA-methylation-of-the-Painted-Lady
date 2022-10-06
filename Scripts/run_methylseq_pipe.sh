#!/bin/bash -l

#Prerequisites: install methylseq pipeline in current directory 
export NXF_HOME="Path to methylseq pipeline directory"

NXF_OPTS='-Xms1g -Xmx4g'


#Example run code, in this cased for the lambda phage. 
nextflow run nf-core/methylseq -name "lambda" --input 'samples/Sample*L*/*R{1,2}*.fastq.gz' -profile uppmax --project "" --max_cpus 20 --max_memory 128GB --fasta lambda_phage.fa --outdir ./van_Lambda --save_reference --comprehensive --cytosine_report --epignome
