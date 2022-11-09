DMR_odds_data <- read.table(file=file.choose(), header=T)
DMR_odds_data$Comparison <- as.factor(DMR_odds_data$Comparison)

library(scales)

DMR_odds_data<-DMR_odds_data[DMR_odds_data$Annotation != "Genes" & DMR_odds_data$Annotation != "Exons",]  

DMR_odds_data$Pvalue_resample <- ifelse(DMR_odds_data$Pvalue_resample <= 0.5, DMR_odds_data$Pvalue_resample*2, (1-DMR_odds_data$Pvalue_resample)*2 )
DMR_odds_data$Significance <- ifelse(DMR_odds_data$Pvalue_resample > 1-(0.025/9) | DMR_odds_data$Pvalue_resample < 0.025/9, "Y", "N" )

DMR_odds_data$Annotation <- gsub("three_prime_UTRs", "3' UTRs", DMR_odds_data$Annotation)
DMR_odds_data$Annotation <- gsub("five_prime_UTRs", "5' UTRs", DMR_odds_data$Annotation)
DMR_odds_data$Annotation <- gsub("DNA_transposons", "DNA transposons", DMR_odds_data$Annotation)

ggplot(DMR_odds_data, aes(x = log10(Odds_ratio), y = reorder(Annotation, Odds_ratio), fill=Comparison,  alpha=Significance)) + 
  scale_alpha_discrete(breaks=c("N", "Y"), range=c(0.3,1) )+
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  scale_fill_manual(name="Comparison", values = c("purple", "lightgreen"))+
  geom_point(size = 3.5, position=position_dodge2(width=0.5), shape=21) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0, 5.5, 1)), labels = seq(0, 5.5, 1),
  #                   limits = log10(c(0.1,5.5))) +
  scale_x_continuous(breaks = log10(seq(0, 10, 1)), labels = seq(0, 10, 1),
                     limits = log10(c(0.1,10)), guide = guide_axis(check.overlap = TRUE)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Odds ratio")+
  theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

