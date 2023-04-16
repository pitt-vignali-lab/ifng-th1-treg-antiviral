# ----------------------------------------------------------------------------
# Revision | Batch correction with simspec and re-clustering
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.1.batchEffCheckAndDataInteg.R
# Input    : 1.Data/revision/1.initDat.RDS
# Output   : 1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS
#            3.Output/revision/2.Integration/
# Upstream : 30_rev_data_init.R
# Note     : cluster_sim_spectrum with label_tag = batch, spectrum_type = corr_ztransform, corr_method = spearman
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/2.Integration")){
  dir.create("3.Output/revision/2.Integration")
}


library(Seurat)
library(tidyverse)
library(simspec)
library(presto)

dat1 <- readRDS("1.Data/treg/TregDat.RDS")
DimPlot(dat1)

# change the smaple names and make them consistent with new settings:
table(dat1$reps)
table(dat1$samples)

# recall the original settings of dat1 from the draft

modDf <- dat1@meta.data[c("RNA_snn_res.0.1", "RNA_snn_res.0.4")]

modDf$RNA_snn_res.0.1 <- as.numeric(as.character(modDf$RNA_snn_res.0.1))

modDf$RNA_snn_res.0.4 <- as.numeric(as.character(modDf$RNA_snn_res.0.4))

modDf$finClu <- as.factor(ifelse(as.numeric(modDf$RNA_snn_res.0.4) == 7, 
                                 4, as.numeric(modDf$RNA_snn_res.0.1)))

dat1 <- AddMetaData(dat1, modDf[-c(1,2)])

Idents(dat1) <- "finClu"
DimPlot(dat1)

annoId <- data.frame(id = 0:4,
                     anno = c("Effector", # 0
                              "Central", # 1
                              "Metabolically \nhyperactive", # 2
                              "Proliferating", # 3
                              "Gene transcription and \nprotein stability" # 4
                     )) |> merge(dat1@meta.data["finClu"] |> 
                                   rownames_to_column(var = "cell"), 
                                 by.x = "id", by.y = "finClu") |> 
  select("cell", "anno") |> 
  column_to_rownames(var = "cell")

dat1 <- AddMetaData(dat1, annoId)

dat1$batch <- "dat1"

DimPlot(dat1, group.by = "anno", label = T)

# dat2 read in 
dat2 <- readRDS("1.Data/revision/1.initDat.RDS")

table(dat2$reps)

dat2$batch <- "dat2"

# merge two datasets
dat <- merge(dat1, y = dat2, add.cell.ids = c("set1", "set2"))

# batch effects checking ----

dat <- NormalizeData(dat)

dat <- FindVariableFeatures(dat)

dat <- ScaleData(dat)

dat <- RunPCA(dat)

dat <- FindNeighbors(dat, dims = 1:20)

dat <- FindClusters(dat, resolution = 0.5)

dat <- RunUMAP(dat, dims = 1:20)

DimPlot(dat, group.by = "batch") +
  coord_fixed()

ggsave("3.Output/revision/2.Integration/batchEffChk.png", 
       width = 6, height = 4, units = "in", dpi = 300)

# remove batch effects ----

dat <- cluster_sim_spectrum(object = dat, 
                            label_tag = "batch", 
                            spectrum_type = "corr_ztransform",
                            corr_method = "spearman")

dat <- RunUMAP(dat, 
               reduction = "css", 
               dims = 1:ncol(Embeddings(dat, "css")))

dat <- FindNeighbors(dat, 
                     reduction = "css", 
                     dims = 1:ncol(Embeddings(dat, "css")))

dat <- FindClusters(dat, 
                    resolution = 0.5)

p1 = DimPlot(dat, label = T) + coord_fixed() + NoLegend()

p2 = DimPlot(subset(dat, batch == "dat1"), group.by = "anno") + coord_fixed()

p1 + p2

ggsave("3.Output/revision/2.Integration/umapPostBatEffRemoval.png", 
       width = 10, height = 4, units = "in")

newAnno <- data.frame(seurat_clusters = c(0, 4, 5, 
                                          1, 2, 3, 7, 
                                          9,
                                          6,
                                          8),
                      newAnno = c(rep("Central", 3),
                                  rep("Effector", 4),
                                  "Gene transcription and \nprotein stability",
                                  "Metabolically \nhyperactive",
                                  "Proliferating"))

dat <- AddMetaData(dat, merge(dat@meta.data["seurat_clusters"] |> 
                                rownames_to_column(var = "cells"),
                              newAnno, by = "seurat_clusters") |> 
                     select("newAnno", "cells") |> 
                     column_to_rownames(var = "cells"))

DimPlot(dat, label = T, group.by = "newAnno") + coord_fixed() + NoLegend()

ggsave("3.Output/revision/2.Integration/umapIntNewDefined.png",
       width = 6, height = 4, units = "in")

DimPlot(dat, group.by = "newAnno", split.by = "samples", ncol = 2) + 
  coord_fixed() 


data.frame(prop.table(table(dat$newAnno, dat$reps), margin = 2)) |> 
  ggplot(aes(x = Freq, y = Var2, fill = Var1)) + 
  geom_bar(stat = "identity")

data.frame(prop.table(table(dat$newAnno, dat$samples), margin = 2)) |> 
  ggplot(aes(x = Freq, y = Var2, fill = Var1)) + 
  geom_bar(stat = "identity")

mkG <- wilcoxauc(dat, group_by = "newAnno") |> 
  filter(padj < 0.05)

write.csv(mkG, "3.Output/revision/2.Integration/markerGenesByCluster.csv")

table(subset(dat, batch == "dat1")$anno, subset(dat, batch == "dat1")$newAnno)

saveRDS(dat, "1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")
