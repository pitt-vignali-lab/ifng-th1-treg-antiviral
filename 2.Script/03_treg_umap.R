# ----------------------------------------------------------------------------
# Setup | Treg UMAP embedding and cluster overview
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.2.formalUMAP.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/2.formalUMAP1*.png
# Upstream : 02_treg_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)

dat <- readRDS("1.Data/treg/TregDat.RDS")

modDf <- dat@meta.data[c("RNA_snn_res.0.1", "RNA_snn_res.0.4")]

modDf$RNA_snn_res.0.1 <- as.numeric(as.character(modDf$RNA_snn_res.0.1))

modDf$RNA_snn_res.0.4 <- as.numeric(as.character(modDf$RNA_snn_res.0.4))

modDf$finClu <- as.factor(ifelse(as.numeric(modDf$RNA_snn_res.0.4) == 7, 
                                 4, as.numeric(modDf$RNA_snn_res.0.1)))

dat <- AddMetaData(dat, modDf[-c(1,2)])

Idents(dat) <- "finClu"

clusterId <- c("Effector", # 0
               "Central", # 1
               "Metabolically\n hyperactive", # 2
               "Proliferating", # 3
               "Gene transcription and\n protein stability" # 4
)

names(clusterId) <- levels(dat)

datLab <- RenameIdents(dat, clusterId)

DimPlot(datLab, label = T) + NoLegend()

df <- data.frame(x = datLab[["umap"]]@cell.embeddings[,1],
                 y = datLab[["umap"]]@cell.embeddings[,2],
                 spl = datLab$samples,
                 lab = datLab@active.ident)

dfLab <- df %>% group_by(lab) %>% 
  summarise(x = mean(x), y = mean(y))

png("3.Output/treg/2.formalUMAP1.png", width = 7, height = 7, units = "in", res = 300)
ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6) +
  geom_label(data = dfLab, aes(x , y, label = lab), col = "black", alpha = 0.1, size = 5.5) +
  theme_classic()+ 
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(legend.position="none",
        text = element_text(size = 15))
dev.off()

df$infe = "Cl. 13 d6"
df[grep("U", df$spl),]$infe = "Na�ve" 

df$spl <- factor(df$spl, levels = c("UC", "UG", "C", "G"))
df$infe <- factor(df$infe, levels = c("Na�ve", "Cl. 13 d6"))

png("3.Output/treg/2.formalUMAP1_split.png", width = 7, height = 7, units = "in", res = 300)
ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6) +
  # geom_label(data = dfLab, aes(x , y, label = lab), col = "black", alpha = 0.1, size = 5.5) +
  theme_classic()+ 
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  facet_wrap(~ infe + spl, nrow = 2, ncol = 2) +
  theme(legend.position="none",
        text = element_text(size = 15),
        strip.text = element_text(size=0,colour="white"),
        strip.background = element_rect(colour="white", fill="white"))
dev.off()

png("3.Output/treg/2.cellTypeBar.png", width = 10, height = 6, units = "in", res = 300)
data.frame(prop.table(table(Idents(datLab), datLab$reps),
                      margin = 2)) %>% 
  ggplot(aes(x = Freq * 100, y = Var2, 
             fill = paste0(Var1, ifelse(Var1 == "Gene transcription and\n protein stability", "   *", "   ns")))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  ylab(NULL) +
  xlab("Percent (%)") +
  theme_bw()+
  scale_fill_discrete(name = "Cell type") +
  theme(legend.position="bottom",
        text = element_text(size=18))
dev.off()

# ----
library(ggsignif)
library(ggpubr)
library(scales)

data <- data.frame(prop.table(table(Idents(dat), dat$reps),
                              margin = 2))

data$samp <- gsub("[0-9]", "", data$Var2)

data$lcmv <- "infected"
data$lcmv[grep("U", data$samp)] <- "uninfected"

data$ifngr1 <- "indeficient"
data$ifngr1[grep("G", data$samp)] <- "deficient"
# data$lcmv <- factor(data$lcmv, levels = c("C", "G"))

# datasub <- data[-grep("U", data$Var2),]


ggplot(data, aes(x = lcmv, y = Freq*100)) + 
  geom_boxplot(aes(fill = Var1), alpha = 0.5) +
  geom_point(aes(col = Var1),position=position_jitterdodge(jitter.width = 0.8),
             size = 3, alpha = 0.2)+
  facet_grid(~Var1)+
  stat_compare_means(aes(group = lcmv),
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
  ggtitle("Effects of LCMV infection") +
  NoLegend()


ggplot(data, aes(x = ifngr1, y = Freq*100)) + 
  geom_boxplot(aes(fill = Var1), alpha = 0.5) +
  geom_point(aes(col = Var1),position=position_jitterdodge(jitter.width = 0.8),
             size = 3, alpha = 0.2)+
  facet_grid(~Var1)+
  stat_compare_means(aes(group = ifngr1),
                     method = "anova",
                     label =  "p.signif",
                     label.x = 1.5)+
  # stat_signif(comparisons=list(c("indeficient", "deficient")),map_signif_level = T, 
  #             test = wilcox.test, textsize = 2, tip_length = 0.005, col = "gray50")+
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
  ggtitle("Effects of Ifngr1 deficiency") + NoLegend()
