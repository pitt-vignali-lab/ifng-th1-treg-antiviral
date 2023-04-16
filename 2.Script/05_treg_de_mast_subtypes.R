# ----------------------------------------------------------------------------
# Analysis | Differential expression within Treg subtypes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/3.2.deMastForSubtypes.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/DEgenes_*.png
# Upstream : 04_treg_de_mast.R
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

DimPlot(dat,
        cells.highlight= WhichCells(dat, expression = Gata3 > 0.1), 
        cols.highlight = "darkred", cols= "grey") + NoLegend()

DimPlot(dat,
        cells.highlight= WhichCells(dat, expression = Tbx21 > 0.1), 
        cols.highlight = "darkblue", cols= "grey") + NoLegend()

datSub <- subset(dat, cell = WhichCells(dat, expression = Tbx21 > 0.1)) 
Idents(datSub) <- "samples"

FindMarkers(datSub, ident.1 = "G", ident.2 = "C",
            test.use = "MAST", logfc.threshold =  0.05) %>% 
  filter(p_val < 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  write.csv("3.Output/treg/deMAST_Tbx21+GvsC.csv")

