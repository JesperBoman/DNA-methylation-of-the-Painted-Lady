#!/usr/bin/env Rscript

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to perform DESeq2 differential expression analysis from the Nextflow rnaseq output
# ==========================================================================================================
# Jesper Boman                      11 nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #



args = commandArgs(trailingOnly=T)

args



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version="3.12")
BiocManager::install("DESeq2")

library(DESeq2)

load("deseq2.dds.RData")

#~Family+Treatment 
dds<-DESeqDataSet(dds, design = ~Group4+Group1)
dds <- DESeq(dds)

expr.res.HDALvLDAL <- results(dds, contrast=c("Group1","HDAL","LDAL"))

expr.res.HDALvLDAL.df <- as.data.frame(expr.res.HDALvLDAL[order(expr.res.HDALvLDAL$pvalue),])
colnames(expr.res.HDALvLDAL.df) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

expr.res.HDALvLDAL.df <- expr.res.HDALvLDAL.df[!is.na(expr.res.HDALvLDAL.df$padj),]
str(expr.res.HDALvLDAL.df[expr.res.HDALvLDAL.df$padj  < 0.1,])


expr.res.HDALvHDLI <- results(dds, contrast=c("Group1","HDAL","HDLI"))

expr.res.HDALvHDLI.df <- as.data.frame(expr.res.HDALvHDLI[order(expr.res.HDALvHDLI$pvalue),])
colnames(expr.res.HDALvHDLI.df) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

expr.res.HDALvHDLI.df <- expr.res.HDALvHDLI.df[!is.na(expr.res.HDALvHDLI.df$padj),]
str(expr.res.HDALvHDLI.df[expr.res.HDALvHDLI.df$padj  < 0.1,])




HDALvLDAL_DMR_genes <- read.table("HDALvLDAL.Genes.geneIDs.list")
colnames(HDALvLDAL_DMR_genes) <- c("gene", "direction")

HDALvHDLI_DMR_genes <- read.table("HDALvHDLI.Genes.geneIDs.list")
colnames(HDALvHDLI_DMR_genes) <- c("gene", "direction")


table(rownames(expr.res.HDALvLDAL.df[expr.res.HDALvLDAL.df$padj  < 0.1,]) %in% HDALvLDAL_DMR_genes$gene) 
HDALvLDAL_DMR_genes[HDALvLDAL_DMR_genes$gene %in% rownames(expr.res.HDALvLDAL.df[expr.res.HDALvLDAL.df$padj  < 0.1,]),]

signDE.HDALvLDAL <- expr.res.HDALvLDAL.df[expr.res.HDALvLDAL.df$padj  < 0.1,]
signDE.HDALvLDAL[rownames(signDE.HDALvLDAL) %in% HDALvLDAL_DMR_genes$gene,]


table(rownames(expr.res.HDALvHDLI.df[expr.res.HDALvHDLI.df$padj  < 0.1,]) %in% HDALvHDLI_DMR_genes$gene) 
HDALvHDLI_DMR_genes[HDALvHDLI_DMR_genes$gene %in% rownames(expr.res.HDALvHDLI.df[expr.res.HDALvHDLI.df$padj  < 0.1,]),]

signDE.HDALvHDLI <- expr.res.HDALvHDLI.df[expr.res.HDALvHDLI.df$padj  < 0.1,]
signDE.HDALvHDLI[rownames(signDE.HDALvHDLI) %in% HDALvHDLI_DMR_genes$gene,]

signDE.HDALvHDLI[rownames(signDE.HDALvHDLI) %in% rownames(signDE.HDALvLDAL),]
signDE.HDALvLDAL[rownames(signDE.HDALvLDAL) %in% rownames(signDE.HDALvHDLI),]

mat <- matrix(c(3, 35, 2091, 13161), nrow=2, ncol=2)
mat <- matrix(c(4, 35, 2021, 13161), nrow=2, ncol=2)
fisher.test(mat)

