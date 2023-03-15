library(ggplot2)
install.packages("ggsignif")
library(ggsignif)

hyper <- c(1606, 1195)
hypo <- c(2965-1606, 2940-1195)

dat <- as.data.frame(cbind(DMRs=c(hyper, hypo), type=c("HDAL higher", "HDAL higher", "HDAL lower", "HDAL lower"), comp=c("HDAL-LDAL", "HDAL-HDLI", "HDAL-LDAL", "HDAL-HDLI")))
dat$DMRs <- as.numeric(dat$DMRs)

df<-data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),y=c(1900, 1900), annotation=c("***", "***"))
dat
positions <- c("HDAL lower", "HDAL higher", "HDAL higher", "HDAL lower")

c("#ef8b41", "#30a2c4", "#30a2c4", "#ef8b41")
ggplot(dat, aes(x=comp, y=DMRs, fill=type))+geom_bar(stat="identity", position="dodge", col="black",width = 0.5)+
  scale_y_continuous(limits = c(0,1900), expand = c(0, 0)) +
  scale_fill_manual(name="Difference", limits=positions, values = c("#30a2c4", "#ef8b41", "#30a2c4", "#ef8b41"))+
  theme_classic()+
  scale_x_discrete(name ="Comparison", labels=c("                                   HDLI <—————> HDAL <—————> LDAL", ""))+
  xlab("Comparison")+
  ylab("Number of DMRs")+
  geom_signif(y_position = c(1800,1800), xmin = c(0.8,1.8), 
              xmax = c(1.2,2.2), annotation = c("***","***"),
              tip_length = 0.03)+
  theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))
