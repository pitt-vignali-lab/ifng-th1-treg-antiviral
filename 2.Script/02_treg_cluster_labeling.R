# ----------------------------------------------------------------------------
# Setup | Assign Treg cluster identities from marker genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.1.cellLabeling.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/
# Upstream : 01_treg_data_init.R
# ----------------------------------------------------------------------------

library(Seurat)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(cluster)
library(clusterProfiler)
library(org.Mm.eg.db)
library(SPIA)
library("EnrichmentBrowser")
library("pathview")
library("org.Mm.eg.db")
library(clustree)

dat <- readRDS("1.Data/treg/TregDat.RDS")

modDf <- dat@meta.data[c("RNA_snn_res.0.1", "RNA_snn_res.0.4")]

modDf$RNA_snn_res.0.1 <- as.numeric(as.character(modDf$RNA_snn_res.0.1))

modDf$RNA_snn_res.0.4 <- as.numeric(as.character(modDf$RNA_snn_res.0.4))

modDf$finClu <- as.factor(ifelse(as.numeric(modDf$RNA_snn_res.0.4) == 7, 
                       4, as.numeric(modDf$RNA_snn_res.0.1)))

dat <- AddMetaData(dat, modDf[-c(1,2)])

Idents(dat) <- "finClu"

DimPlot(dat)

DimPlot(dat, split.by = "samples")

markers <- FindAllMarkers(dat, only.pos = T)
markers_sig <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

write.csv(markers_sig, "3.Output/treg/top20SigMarker.csv")

# Idents(dat) <- "RNA_snn_res.0.1"

# naive CCR7hiTCF7hiHLA-DRlowFOXP3low 1
# 
# cTreg 
VlnPlot(dat, features = c("Tcf7", "Bcl2", "Sell", "S1pr1"), ncol = 4)

# eTreg 
VlnPlot(dat, features = c("Tnfrsf9", "Relb", "Ikzf2", "Pdcd1"), ncol = 4)

# cell cycle 3
mk <- FindAllMarkers(dat, only.pos = T)

ids=bitr(mk$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')

sce.markers=merge(mk,ids,by.x='gene',by.y='SYMBOL')

sce.markers <- sce.markers[sce.markers$p_val_adj < 0.05, ]
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="mmu", pvalueCutoff=0.05)

pDf <- xx@compareClusterResult

ggplot(pDf[# pDf$ID %in% names(table(pDf$ID)[table(pDf$ID) == 2])&
             pDf$Cluster %in% c(2, 3, 4),],
              aes(x = Cluster, y = Description)) + 
  geom_point()
  


