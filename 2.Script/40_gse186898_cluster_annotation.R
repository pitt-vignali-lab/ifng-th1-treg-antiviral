# ----------------------------------------------------------------------------
# Data mining | Annotate immune cell types in GSE186898 from marker genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/6.3.2.gse186898ClustAnno.R
# Input    : 1.Data/revision/
# Output   : 3.Output/revision/
# Upstream : 39_gse186898_data_init.R
# Note     : Cell types assigned with the markers reported in the source publication
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/")){
  dir.create("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/")
}

library(Seurat)
library(tidyverse)
library(presto)

dat <- readRDS("3.Output/revision/6.Chemokins/GSE186898/1.initDat.RDS")

dat <- FindNeighbors(dat, dims = 1:14)
dat <- RunUMAP(dat, dims = 1:14)


dat <- FindClusters(dat, resolution = 0.15)

DimPlot(dat, reduction = 'umap', label = T) +
  coord_fixed()
ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/originalClusterDistr.png")

DotPlot(dat, features = c("Ighd", "Trbc2", "Ace", "C1qb", "Fscn1",
                          "S100a9", "Stmn1", "Jchain", "Klrb1c", "Mcpt8",
                          "Foxp3", "Cd8a", "Cd4")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/originalMarkExp.png",
       bg = "white")

datSub <- subset(dat, idents = 3)
datSub <- FindNeighbors(datSub, dims = 1:14)
datSub <- RunUMAP(datSub, dims = 1:14)
datSub <- FindClusters(datSub, resolution = 0.15)

DimPlot(datSub, label = T) + 
  coord_fixed()
ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/subclstCD4UMAPGetTreg.png")

DotPlot(datSub, features = c("Foxp3")) +
  coord_fixed()
ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/foxp3ExpInSubclst.png",
       bg = "white")

DimPlot(dat, cells.highlight = Cells(subset(datSub, idents = 2)))+
  coord_fixed()
ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/tregDistrInAll.png")


marks <- wilcoxauc(dat)
write.csv(marks, "3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/markInOrigClst.csv")

marks %>% filter(group == 10) %>% arrange(padj) %>% 
  head

annoDf <- data.frame(idents = 0:10,
                     anno = c("B cell", "Monocyte", "Macrophage",
                              "CD4 T cell", "CD8 T cell", "Dendritic",
                              "Neutrophil", "B cell", "Plasma", 
                              "Neutrophil", "B cell")) %>% 
  merge(dat@meta.data['seurat_clusters'] %>% rownames_to_column("cells"),
        by.x = "idents", by.y = "seurat_clusters")

annoDf[annoDf$cells %in% Cells(subset(datSub, idents = 2)),]$anno <- "Treg"

dat <- AddMetaData(dat, annoDf %>% select("anno", "cells") %>%
                     column_to_rownames("cells"))

DimPlot(dat, group.by = "anno", label = T, label.size = 5) +
  coord_fixed() +
  NoLegend()

ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/annoClstDistr.png",
       width = 6, height = 8)

DotPlot(dat, features = c("Ighd", "Trbc2", "Ace", "C1qb", "Fscn1",
                          "S100a9", "Stmn1", "Jchain", "Klrb1c", "Mcpt8",
                          "Foxp3", "Cd8a", "Cd4"),
        group.by = "anno") +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/annoClstMarks.png",
       width = 8, height = 6, bg = "white")

wilcoxauc(dat, group_by = "anno") %>% 
  write.csv("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/annoClstMarksList.csv")

Embeddings(dat, reduction = "umap") %>% 
  merge(dat@meta.data[c("title", "anno")], by = 0) %>% 
  separate(title, into = c("geno", "time"), sep = "_") %>% 
  mutate(time = factor(time, levels = c("d1", "d9", "d25", "d30"))) %>% 
  ggplot(aes(UMAP_1, UMAP_2, col = anno)) +
  geom_point() +
  facet_grid(geno ~ time) +
  theme_bw() +
  coord_fixed()
ggsave("3.Output/revision/6.Chemokins/GSE186898/2.ClusterAnnotation/annoClstSplited.png",
       width = 10, height = 6, bg = "white")

saveRDS(dat, "3.Output/revision/6.Chemokins/GSE186898/2.datAnnotated.RDS")
