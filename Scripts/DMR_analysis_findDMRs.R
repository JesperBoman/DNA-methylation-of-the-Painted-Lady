#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2022-03-15
#### Following tutorials that can be find under vignette: https://www.bioconductor.org/packages/release/bioc/html/bsseq.html




library("bsseq")
library("BiocParallel")

args = commandArgs(trailingOnly = TRUE)

load(args[1])
load(args[2])


BS.cov <- getCoverage(BS_data.fit)

#Filter for loci in which at least 2 samples have at least 2 reads in both sample groups
#Assumes two groups with sample size of 3, change accordingly if your sample size is different
keepLoci<- which(rowSums(BS.cov[, 1:3] >= 2) >= 2 &
                            rowSums(BS.cov[, 4:6] >= 2) >= 2)
length(keepLoci)

comp.fit.filt <- BS_data.fit[keepLoci,]


comp.fit.filt.tstat <- BSmooth.tstat(comp.fit.filt, 
                                          group1 = colnames(BS_data)[1:3],
                                          group2 = colnames(BS_data)[4:6], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "t_plot.png")
plot(comp.fit.filt.tstat )
dev.off()


dmrs0 <- dmrFinder(comp.fit.filt.tstat , qcutoff = c(0.01, 0.99))

#Filter for dmrs in which at least 3 CpG loci show an absolute methylation difference of at least 0.1
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(comp.fit.filt, file="comp.fit.filt.rda")
save(dmrs, file="dmrs.rda")
write.table(dmrs, file="dmrs.txt", sep="\t", quote =F, row.names=F)
dmrs$start <- dmrs$start - 1 
write.table(dmrs, file="dmrs.bed", sep="\t", quote =F, row.names=F)
