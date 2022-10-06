#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis.R *input_files*
#### Jesper Boman - 2020-09-05

args = commandArgs(trailingOnly = TRUE)

library("bsseq")
library("BiocParallel")

#Reading data and creating BSseq object

for (i in 1:length(args)){
BS_table <- read.table(file=args[i])

sample_name_temp <- unlist( strsplit(args[i], "\\.") )[1]
        if(i == 1){
                BS_data <-BSseq(chr=BS_table$V1, pos=BS_table$V2, M=matrix(BS_table$V3), Cov=matrix(BS_table$V4), sampleNames=sample_name_temp)
                }
        else{
                BS_data <- combine(BS_data, BSseq(chr=BS_table$V1, pos=BS_table$V2, M=matrix(BS_table$V3), Cov=matrix(BS_table$V4), sampleNames=sample_name_temp))
                }
}
save(BS_data, file="BS_data.rda")


#Smoothing

#Set the number of cores here
options(mc.cores=14)

BS_data.fit <- BSmooth(BSseq = sort(BS_data), BPPARAM = MulticoreParam(workers = 14, progressbar = TRUE))

save(BS_data.fit, file="BS_data.fit.rda")
