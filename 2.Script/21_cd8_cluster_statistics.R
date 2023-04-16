# ----------------------------------------------------------------------------
# ED Fig 9b-c | Cluster composition statistics by genotype
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.2.cellClusterStatistic.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggpubr)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")

Idents(dat) <- "RNA_snn_res.0.3"

DimPlot(dat)

data <- data.frame(prop.table(table(Idents(dat), dat$reps),
                              margin = 2))

data$samp <- gsub("[0-9]", "", data$Var2)

data$ifngr1 <- "indeficient"
data$ifngr1[grep("G", data$samp)] <- "deficient"
# data$lcmv <- factor(data$lcmv, levels = c("C", "G"))

# datasub <- data[-grep("U", data$Var2),]


data.frame(prop.table(table(Idents(dat), dat$reps),
                      margin = 2)) %>% 
  rename(CellType = Var1, Sample = Var2) %>% 
  ggplot(aes(y = Sample, x = Freq*100, fill=CellType, label = paste0(round(Freq, 2) * 100,"%"))) +
  geom_bar(stat='identity', colour="black", alpha = 0.5) +
  xlab("Persentage (%)") +
  theme_bw()


ggplot(data, aes(x = ifngr1, y = Freq*100)) + 
  geom_boxplot(aes(fill = Var1), alpha = 0.5) +
  geom_point(aes(col = Var1),position=position_jitterdodge(jitter.width = 0.8),
             size = 3, alpha = 0.2)+
  facet_grid(~Var1)+
  stat_compare_means(aes(group = ifngr1),
                     method = "anova",
                     label =  "p.signif",
                     label.x = 1.5) +
  # stat_signif(comparisons=list(c("uninfected", "infected")),map_signif_level = T,
  #             test =  wilcox.test, textsize = 2, tip_length = 0.005, col = "gray50")+
  scale_y_log10()+
  ylab("Percentage (%)")+
  xlab(NULL)+
  # stat_compare_means(comparisons=list(c('C','G'), c('C', 'UC'), c('G', "UG")))+
  # stat_compare_means(label.y = 0.8)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        # strip.text.x = element_blank(),
        strip.background.x = element_blank()
  )+ 
  ggtitle("Effects of Ifngr1 deficiency") +
  NoLegend()
