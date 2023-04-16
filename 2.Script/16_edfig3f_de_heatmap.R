# ----------------------------------------------------------------------------
# ED Fig 3f | Heat map of significantly differentially expressed genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/S2C(heatmap).R
# Input    : 1.Data/treg/TregDat.RDS
#            3.Output/treg/deGsForVolcanoPlot.csv
# Output   : 3.Output/treg/S2CDenSel.png
# Upstream : 04_treg_de_mast.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dat <- readRDS("1.Data/treg/TregDat.RDS")

# Idents(dat) <- "samples"

# old ----
# deGs <- FindMarkers(dat, ident.1 = "C", ident.2 = "G")
# deGs <- deGs[deGs$p_val_adj < 0.05,]

# MAST ----

deGs <- deOverall <- read.csv("3.Output/treg/deGsForVolcanoPlot.csv")

# heat maps ----
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(extrafont)

loadfonts(device = "win", quiet = TRUE)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

dfhm <- AverageExpression(dat,features = deGs[deGs$p_val_adj<0.05,]$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

htDeGs <- c(head(deGs[deGs$avg_log2FC > 0 & sort(deGs$p_val, decreasing = F),]$gene, 20),
            head(deGs[deGs$avg_log2FC < 0 & sort(deGs$p_val, decreasing = F),]$gene, 20))

htGenes <- c('Ifngr1', 'Ly6c2', 'Cxcr3', 'AW112010', 'Stat1', 'Tbx21', 
             'Il7r', 'Gata3', 'Klrg1', 'Nkg7', 'Ifit3')

htGenes <- unique(c(htGenes,
                    htDeGs))

ha = rowAnnotation(foo = anno_mark(at = which(rownames(htmp) %in%  htGenes), 
                                   labels_gp = gpar(fontsize = 18, fontface = 2),
                                   labels = rownames(htmp)[rownames(htmp) %in%  htGenes]))

png("3.Output/treg/S2CDenSel.png", height = 12.5, width = 10, units = "in", res = 300)
draw(Heatmap(htmp,
        na_col = "black",
        rect_gp = gpar(col = "gray50", lwd = 0.8),
        row_dend_width = unit(2, "cm"),
        col = col_fun,
        right_annotation = ha, 
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = T,
        show_row_dend = T,
        cluster_columns = F,
        column_split = c(1, rep(2, 3), 3, rep(4, 3)),
        # column_names_side = "bottom",
        # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
        # column_names_rot = 0,
        # column_names_gp = gpar(fontsize = 10),
        column_title = NULL,
        row_names_side = "left",
        show_column_names = F,
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
                                                                             fontsize = 10))),
        heatmap_legend_param = list(direction = "horizontal",
                                    title_position = "topcenter",
                                    legend_width = unit(7.5, "cm"))), heatmap_legend_side = "bottom")

dev.off()

