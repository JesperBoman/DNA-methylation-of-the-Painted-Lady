#!/usr/bin/env Rscript

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==========================================================================================================
# A script to perform gene ontology analysis
# Based on original script by Karin Näsvall
# https://github.com/DaSh-bash/GenomeVanessa/blob/main/Badirate_incl_topGO_220114/topGO/scripts/topGO_Vanessa_cardui_gains.Rmd
# ==========================================================================================================
# Jesper Boman                      9 nov 2022
# ==========================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #



BiocManager::install("topGO")
library(topGO)
library(ggplot2)

geneID2GO <- readMappings(file="vanCar_GO_genes.list", sep = "\t", IDsep = ",")

str(geneID2GO)
geneNames <- names(geneID2GO)



#read the list of genes of interest
candidate_genes <- read.table("HDALvHDLI.CDS.geneIDs.list", header = F)
#candidate_genes <- rownames(expr.res.HDALvHDLI.df[expr.res.HDALvHDLI.df$padj  < 0.1,])
colnames(candidate_genes) <- c("geneid", "direction")

head(candidate_genes)
candidate_genes <- as.character(candidate_genes$geneid)
geneList <- factor(as.integer(geneNames %in% candidate_genes))
names(geneList) <- geneNames
head(geneList)
str(geneList)
#create the topGOdata object for each ontology
GO_data_BP <- new("topGOdata", 
                  ontology="BP", 
                  allGenes=geneList, 
                  annot=annFUN.gene2GO, 
                  gene2GO=geneID2GO, 
                  nodeSize=1)
GO_data_BP


#"when only a list of interesting genes is provided, 
#the user can use only tests statistics that are based on gene counts, 
#like Fisher’s exact test, Z score and alike."

##Algorithms:
#classic uses all GO terms separately
#elim eliminates genes from parentGOterms if child is more specific (bottom up analysis)
#weight trying to determin if GOterm is better representing the list of interesting genes (is more enriched) than any other node from its neighborhood
#parentChild 

#biological process
BP_resultFisher_weight01 <- runTest(GO_data_BP, algorithm = "weight01", statistic = "fisher")
BP_resultFisher_weight01

BP_resultFisher_classic <- runTest(GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_resultFisher_classic

BP_resultFisher_parentChild <- runTest(GO_data_BP, algorithm = "parentChild", statistic = "fisher")
BP_resultFisher_parentChild

BP_resultFisher_elim <- runTest(GO_data_BP, algorithm = "elim", statistic = "fisher")
BP_resultFisher_elim



#allGO = usedGO(object = GOdata)
#topNodes = length(allGO) in GenTable() to get a table with all GO:s to do fdr
allGO = usedGO(object = GO_data_BP)
allRes_BP <- GenTable(GO_data_BP, 
                      weight01Fisher = BP_resultFisher_weight01, 
                      classicFisher = BP_resultFisher_classic,
                      parentChFisher=BP_resultFisher_parentChild,
                      elimFisher=BP_resultFisher_elim,
                      orderBy = "weight01Fisher", 
                      ranksOf = "weight01Fisher", 
                      topNodes = length(allGO),
                      numChar=1000)
allRes_BP$weight01Fisher <- as.numeric(allRes_BP$weight01Fisher)
allRes_BP[allRes_BP$weight01Fisher<=0.05, ]
allRes_BP$parentChFisher <- as.numeric(allRes_BP$parentChFisher)
#allRes_BP[allRes_BP$parentChFisher<=0.05, ]
allRes_BP$elimFisher <- as.numeric(allRes_BP$elimFisher)
#allRes_BP[allRes_BP$elimFisher<=0.05, ]

#multiple test correction, method Benjamini-Hochberg
allRes_BP$p_adj <- p.adjust(allRes_BP$weight01Fisher, method = "BH")
allRes_BP$p_adj_pc <- p.adjust(allRes_BP$parentChFisher, method = "BH")
allRes_BP$p_adj_elim <- p.adjust(allRes_BP$elimFisher, method = "BH")

#adj p-value below 0.1
allRes_BP[allRes_BP$p_adj<0.1, ]
allRes_BP[allRes_BP$p_adj_pc<0.1, ]
allRes_BP[allRes_BP$p_adj_elim<0.1, ]
allRes_BP$GO_class <- "BP"
