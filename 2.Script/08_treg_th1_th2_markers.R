# ----------------------------------------------------------------------------
# Fig 3 | Feature and violin plots of TH1 and TH2 marker genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/5.1.th1th2Markers.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/Th1HallmarkOverall.png
#            3.Output/treg/Th2HallmarkOverall.png
# Upstream : 02_treg_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggalt)
library(ggforce)
library(ggpubr)
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
DimPlot(datLab)

# Th1
FeaturePlot(dat, features = c("Tbx21", "Ifng", "Il12rb1", "Rela"), label = T)

FeaturePlot(datLab, features = c("Gata3", "Il2ra", "Il4ra", "Il4", "Il2", "Maf",
                                 "Lima1", "Auh", "Runx2", "Smad2"))

Idents(dat) <- "RNA_snn_res.0.5"

DimPlot(dat, label = T)
datSub <- subset(dat, idents = 12)

DimPlot(dat,
        cells.highlight= WhichCells(dat, expression = Ifng > 0.5), 
        cols.highlight = "darkred", cols= "grey")

th1Df <- data_frame()
for (i in c("Tbx21", "Ifng", "Il12rb1", "Rela")) {
  p <- VlnPlot(dat, features = i,
          group.by = "samples")
  p$data$gene <- names(p$data)[1]
  names(p$data)[1] <- "value"
  th1Df <- rbind(th1Df, p$data)
}

th1Df <- th1Df[th1Df$ident %in% c("G", "C"),]
th1Df$ident <- ifelse(th1Df$ident == "G", "Mutant", "Control")

png("3.Output/treg/Th1HallmarkOverall.png", width = 5, height = 4, units = "in", res = 300)
ggplot(th1Df[th1Df$value > 0.1, ], aes(x = ident, y = value, fill = ident)) +
  geom_boxplot() +
  stat_compare_means(aes(group = ident),
                     method = "anova",
                     label =  "p.signif",
                     label.x = 1.5) +
  scale_fill_manual(values = c("#1B9E77","#D95F02")) +
  facet_grid(~gene) +
  ylab("Expression level") +
  xlab(NULL) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

Idents(datSub) <- "samples"


deG <- FindMarkers(datSub, ident.1 = "G", ident.2 = "C", test.use = "MAST", logfc.threshold = "0.05") 

write.csv(deG, "3.Output/treg/deMAST_th1like.csv")

df <- data.frame(x = dat[["umap"]]@cell.embeddings[,1],
                 y = dat[["umap"]]@cell.embeddings[,2],
                 lab = dat@active.ident)

png("3.Output/treg/umapTh1Highlight.png", width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x, y)) +
  geom_point(col = "gray60")+
  geom_point(data = df[df$lab == 12,], aes(x, y), col = "tomato", alpha = 0.2) +
  geom_mark_hull(data = df[df$lab == 12,], aes(x, y), expand=0.01,radius=0.01,
                 linetype="dotted") +
  theme_void()
dev.off()

th2Df <- data_frame()
for (i in c("Gata3", "Il2ra", "Il4ra", "Il4", "Il2", "Maf","Lima1", "Auh", "Runx2", "Smad2")) {
  p <- VlnPlot(dat, features = i,
               group.by = "samples")
  p$data$gene <- names(p$data)[1]
  names(p$data)[1] <- "value"
  th2Df <- rbind(th2Df, p$data)
}

th2Df <- th2Df[th2Df$ident %in% c("G", "C"),]
th2Df$ident <- ifelse(th2Df$ident == "G", "Mutant", "Control")

png("3.Output/treg/Th2HallmarkOverall.png", width = 12.5, height = 5, units = "in", res = 300)
ggplot(th2Df[th2Df$value > 0.1, ], aes(x = ident, y = value, fill = ident)) +
  geom_boxplot() +
  stat_compare_means(aes(group = ident),
                     method = "anova",
                     label =  "p.signif",
                     label.x = 1.5) +
  scale_fill_manual(values = c("#1B9E77","#D95F02")) +
  facet_grid(~gene) +
  ylab("Expression level") +
  xlab(NULL) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()
