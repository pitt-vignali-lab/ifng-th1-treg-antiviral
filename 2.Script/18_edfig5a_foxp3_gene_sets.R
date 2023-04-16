# ----------------------------------------------------------------------------
# ED Fig 5a | Expression of Foxp3-induced and Foxp3-repressed gene sets
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/targetGentIndheatmap.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/foxp3Genesheat.png
# Upstream : 04_treg_de_mast.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(ggpubr)

dat <- readRDS("1.Data/treg/TregDat.RDS")

modDf <- dat@meta.data[c("RNA_snn_res.0.1", "RNA_snn_res.0.4")]

modDf$RNA_snn_res.0.1 <- as.numeric(as.character(modDf$RNA_snn_res.0.1))

modDf$RNA_snn_res.0.4 <- as.numeric(as.character(modDf$RNA_snn_res.0.4))

modDf$finClu <- as.factor(ifelse(as.numeric(modDf$RNA_snn_res.0.4) == 7, 
                                 4, as.numeric(modDf$RNA_snn_res.0.1)))
dat <- AddMetaData(dat, modDf[-c(1,2)])

Idents(dat) <- "finClu"

datSub <- subset(dat, idents = 0)

Idents(dat) <- "samples"

gDe <- FindMarkers(dat, ident.1 = "G", ident.2 = "C")
gDe <- gDe[gDe$p_val_adj < 0.05,]

# Gene read in ----
gIndu <- strsplit("Zscan29
Cd200r1
Nrp1
Foxp3
Pcyt1a
Gpr83
Ccr4
Tgm2
Itgae
Neb
Icos
Gbp2
Lrrc32
Izumo1r
Tnfrsf4
Tnfrsf18
Il2ra
Ctla4
Socs2
P2ry10
Il10ra
Nt5e
Pglyrp1
Ikzf4
Dusp4
Gbp3
Il2rb
Lag3
Capg
Prdm1
Arhgap31
Igf1r
Tgfbr1
Csf1
Il6ra
Jun
Tnf
Lta
Gnaq
Nrn1
Sdc4
Cd80
Maf
Cd69
Nr4a1
Il1rl1
Fos
Rorc
Gm4951", "\n")[[1]]

gRepr <- strsplit("Vipr1
Ikzf2
Slco4a1
Ccr1
Fosl2
Pdcd1
Lgals1
Ldha
Slc16a3
Egr1
Cxcr4
Il4
Il5
Anln
Pole
Ctla2a
Ptpn22
Mettl9
Pde3b
Igfbp4
Cd226
Gen1
Ezh2
Txn1
Tbx21
Id2
Gzma
Eomes
Il2
Fgl2
Ccr8
Cxcr6
Il18rap
Il12rb2
Lgals3
Gzmb
Xcl1", "\n")[[1]]

# heatmap ----
deGs <- data.frame(gene = c(gIndu, gRepr),
                   cat = c(rep("Foxp3 Induced", length(gIndu)),
                           rep("Foxp3 Repressed", length(gRepr))))

deOverall <- read.csv("3.Output/treg/deGsForVolcanoPlot.csv")

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

rownames(htmp) <- paste0(rownames(htmp), 
                         " (adj_p: ", 
                         format(sapply(rownames(htmp), 
                                       function(x) deOverall[deOverall$gene == x,]$p_val_adj), 
                                digits = 3), 
                         ")")

png("3.Output/treg/foxp3Genesheat.png.png", height = 15, width = 11, units = "in", res = 300)
draw(Heatmap(htmp,
             na_col = "black",
             rect_gp = gpar(col = "gray50", lwd = 0.8),
             row_dend_width = unit(2, "cm"),
             col = col_fun,
             row_names_gp = gpar(fontsize = 12, fontface = 2,
                                 col = c("#7C4935", "#457373")),
             cluster_row_slices = FALSE, 
             show_row_dend = T,
             cluster_columns = F,
             column_split = c(1, rep(2, 3), 3, rep(4, 3)),
             # split = deGs$cat,
             row_split = c(rep("Indu", length(gIndu)),
                           rep("Repres", length(gRepr))),
             column_title = NULL,
             row_names_side = "left",
             show_column_names = F,
             row_title = NULL,
             name = "Normalized Gene Expression",
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#fb0707", "#fb0707",
                                                                                    "#7434a4", "#7434a4")),
                                                                 labels = c(
                                                                   expression(Uninfected-""), 
                                                                   expression(Infected-Foxp3^{Cre-YFP}),
                                                                   expression(Uninfected-""),
                                                                   expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                                 ),
                                                                 labels_gp = gpar(col = "white", 
                                                                                  fontsize = 10,
                                                                                  fontface = 2))),
             
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#7C4935", "#457373")),
                                                              labels = c("Foxp3 induced genes", "Foxp3 repressed genes"),
                                                              labels_gp = gpar(col = "white", 
                                                                               fontsize = 10,
                                                                               fontface = 2))),
             
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = unit(7.5, "cm"))), heatmap_legend_side = "bottom")
dev.off()
