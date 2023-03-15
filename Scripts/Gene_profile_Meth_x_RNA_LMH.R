#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_plots_comb.R input_data 
#### Jesper Boman - 2021-05-25

library(ggplot2)
library(plyr)


args = commandArgs(trailingOnly=T)

args="meth_gene_5kb_flanks_100bp_seg_data_coords_expr"

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")


colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Sample", "Group", "Segment_type", "Segment_number", "Coordinate", "TPM")

data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149


data$Expression_cat <- NA
for (ind in unique(data$Sample)){
  
  data[data$TPM <= quantile(data[data$Sample == ind,]$TPM, 1/10, na.rm=T) & data$Sample == ind,]$Expression_cat <- "L"
  data[data$TPM > quantile(data[data$Sample == ind,]$TPM, 1/10, na.rm=T) & data$TPM <= quantile(data[data$Sample == ind,]$TPM, 9/10, na.rm=T)  & data$Sample == ind,]$Expression_cat <- "M"
  data[data$TPM > quantile(data[data$Sample == ind,]$TPM, 9/10, na.rm=T) & data$TPM <= quantile(data[data$Sample == ind,]$TPM, 1, na.rm=T) & data$Sample == ind,]$Expression_cat <- "H"
  
}

data<-data[!is.na(data$Mean_per_dinuc),]

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "Sample", "Expression_cat"),function(x) mean(x$Mean_per_dinuc, na.rm=T))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)

res_body <- ddply(subData, c("Segment_number", "Segment_type", "Sample", "Expression_cat"),function(x) mean(x$Mean_per_dinuc, na.rm=T))


all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "Sample", "Expression_cat", "Methylation_level")
all.res$Segment_number <- as.integer(all.res$Segment_number)
all.res2 <- all.res[!is.na(all.res$Expression_cat),]

all.res2$Expression_cat <- factor(all.res2$Expression_cat, levels=c("L", "M", "H"))

ggplot(data=all.res2, aes(x=Segment_number, y=Methylation_level*100, col=Expression_cat, group=interaction(Sample, Expression_cat))) +
  theme_classic()+
  geom_line() +
  scale_color_manual(name="Expression category", values = c("#1E88E5", "#FFC107", "#D81B60"), labels = c("Low", "Medium", "High")) +
  ylab("Methylation level (%)")+
  ylim(0,18)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="top", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=14),  legend.title=element_text(size=14))

#Local polynomial regression
ggplot(data=all.res2, aes(x=Segment_number, y=Methylation_level*100, col=Expression_cat, fill=Expression_cat, group=interaction(Sample, Expression_cat))) +
  theme_classic()+
  #geom_line() +
  geom_smooth(method = 'loess', formula = 'y ~ x', alpha=0.2, span=0.2)+
  scale_color_manual(name="Expression category", values = c("#1E88E5", "#FFC107", "#D81B60"), labels = c("Low", "Medium", "High")) +
  scale_fill_manual(name="Expression category", values = c("#1E88E5", "#FFC107", "#D81B60"), labels = c("Low", "Medium", "High")) +
  ylab("Methylation level (%)")+
  ylim(0,18)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="top", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=14),  legend.title=element_text(size=14))

