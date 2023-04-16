# ----------------------------------------------------------------------------
# ED Fig 9a | CD8 T cell UMAP embedding
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.3.FormalUMAP.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/2.formalUMAP*.png
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
require(scales)
library(ggpubr)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")
DefaultAssay(dat) <- "RNA"
Idents(dat) <- "RNA_snn_res.0.3"

datSub <- subset(dat, idents = 8, invert = TRUE)

clusterId <- c("Activated", # 0
               "Cytotoxic", # 1
               "Progenitor\n exhausted", # 2
               "Proliferating", # 3
               "Exhausted", # 4
               "Memory", # 5
               "IFN signaling", # 6
               "Metabolically\n hyperactive" # 7
)

names(clusterId) <- levels(datSub)

datLab <- RenameIdents(datSub, clusterId)

DimPlot(datLab, label = T) + NoLegend()


df <- data.frame(x = datLab[["umap"]]@cell.embeddings[,1],
                 y = datLab[["umap"]]@cell.embeddings[,2],
                 lab = datLab@active.ident,
                 rep = datLab$reps)

dfLab <- df %>% group_by(lab) %>% 
  summarise(x = mean(x), y = mean(y))

png("3.Output/cd8/2.formalUMAP2.png", width = 7, height = 7, units = "in", res = 300)
ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6) +
  geom_label(data = dfLab, aes(x , y, label = lab), col = "black", alpha = 0.2, size = 6) +
  theme_classic()+
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(legend.position="none")
dev.off()

png("3.Output/cd8/2.formalUMAP_G1-3.png", width = 15, height = 5, units = "in", res = 300)
ggplot(df[df$rep %in% c("G1", "G2", "G3"),], aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6) +
  # geom_label(data = dfLab, aes(x , y, label = lab), col = "black", alpha = 0.5) +
  facet_grid(~rep)+
  theme_void()+ 
  theme(legend.position="none")

dev.off()

data.frame(prop.table(table(Idents(datLab), datLab$reps),
                      margin = 2)) %>% 
  rename(CellType = Var1, Sample = Var2) %>% 
  ggplot(aes(y = Sample, x = Freq*100, fill=CellType, label = paste0(round(Freq, 2) * 100,"%"))) +
  geom_bar(stat='identity', colour="black", alpha = 0.5) +
  xlab("Persentage (%)") +
  theme_bw()+ 
  theme(legend.position="bottom")

data <- data.frame(prop.table(table(Idents(datLab), datLab$reps),
                              margin = 2))
data$ifngr1 <- "indeficient"
data$ifngr1[grep("G", data$Var2)] <- "deficient"

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

## Bar chart
dataBar <- data[data$Var1 == "Memory",] %>% 
  group_by(ifngr1) %>% 
  summarise(mean = mean(Freq),
            sd = sd(Freq))

ggplot(data[data$Var1 == "Memory",], aes(col = ifngr1, fill = ifngr1))  +
  geom_point(aes(x = ifngr1, y = Freq*100),
             position=position_jitterdodge(jitter.width = 0.8),
             size = 3, alpha = 0.5)+
  stat_compare_means(aes(x = ifngr1, y = Freq*100, group = ifngr1), 
                     method = "anova",
                     label.y = 5.3) +
  stat_compare_means(aes(x = ifngr1, y = Freq*100, group = ifngr1),
                     method = "anova",
                     label =  "p.signif",
                     label.x = 1.5) + 
  geom_bar(data = dataBar, aes(x = ifngr1, y = mean*100), 
           stat = "identity",
           alpha = 0.5)+
  geom_errorbar(data = dataBar, 
                aes(x=ifngr1, ymin=(mean-sd)*100, ymax=(mean+sd)*100), 
                width=0.4, colour="gray30", alpha=0.9, size=0.5) +
  # scale_y_log10()+
  ylab("Percentage (%)")+
  xlab(NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        # strip.text.x = element_blank(),
        strip.background.x = element_blank()
  )+ 
  ggtitle("Effects of Ifngr1 mutation") +
  NoLegend()
