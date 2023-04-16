# ----------------------------------------------------------------------------
# Data mining | Linkage view of the cytokine-receptor network
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/6.3.3.gse186898celltalkerLinkage.R
# Input    : 3.Output/revision/
# Output   : 3.Output/revision/
# Upstream : 41_gse186898_celltalker.R
# ----------------------------------------------------------------------------

library(Seurat)
library(presto)
library(tidyverse)

lrDf <- read.csv("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/Overall_lr.csv")
datTreg <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")
datCd8 <- readRDS("1.Data/cd8/cd8Dat.RDS")

lrTreg <- unique(rbind(data.frame(gene = lrDf[lrDf$cell_type1 == "Treg",]$ligand,
                                  type = "ligand"),
                       data.frame(gene =lrDf[lrDf$cell_type2 == "Treg",]$receptpor,
                                  type = "receptor")))

lrCd8 <- unique(rbind(data.frame(gene = lrDf[lrDf$cell_type1 == "CD8 T cell",]$ligand,
                                 type = "ligand"),
                      data.frame(gene =lrDf[lrDf$cell_type2 == "CD8 T cell",]$receptpor,
                                 type = "receptor")))


deTreg <- wilcoxauc(datTreg, group_by = "samples", groups_use = c("C", "G")) %>% 
  filter(group == "G",
         feature %in% lrTreg$gene)

deCd8 <- wilcoxauc(datCd8, group_by = "samples", groups_use = c("C", "G")) %>% 
  filter(group == "G",
         feature %in% lrCd8$gene)



DimPlot(datTreg, split.by = "samples") + 
  coord_fixed()

DimPlot(datCd8, split.by = "samples") +
  coord_fixed()


