# ----------------------------------------------------------------------------
# Fig 3a | Heat map of TH1 and TH2 gene expression
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/Th1Th2SpectrumChart.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/
# Upstream : 02_treg_cluster_labeling.R
# Note     : ComplexHeatmap
#            batches scaled separately to avoid over-correction
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(ggpubr)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))


dat <- readRDS("1.Data/treg/TregDat.RDS")

Idents(dat) <- "samples"

# deGs <- FindMarkers(dat, ident.1 = "C", ident.2 = "G")
# deGs <- strsplit("Etv5, Hopx, Onecut2, Tbx21/ Ifngr1, Il18r1, Cxcr6, Il18rap, Klrc1, Il12rb2/ Ifng, Il34, Csf1", ", |/ ")[[1]]
# deGs <- strsplit("Batf, Epas1, Gata3, Nfil3/ Fcgrt/ Il4, Tnfsf13b", ", |/ ")[[1]]

# barplot ----
deGs <- strsplit("Tbx21/ Ifngr1, Il18r1, Il18rap, Il12rb2/ Ifng", ", |/ ")[[1]]

deGs <- data.frame(gene = deGs, cat = c(rep("Transcription factor", 1), 
                                        rep("Receptors", 4),
                                        rep("Cytokiunes", 1)))

deGs <- strsplit("Batf, Gata3/ Il4", ", |/ ")[[1]]
deGs <- data.frame(gene = deGs, cat = c(rep("Transcription factor", 2), 
                                        # rep("Receptors", 1),
                                        rep("Cytokiunes", 1)))

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- as.data.frame(dfhm) %>% 
  rownames_to_column("gene") %>% 
  reshape::melt(id = c("gene")) %>% 
  merge(deGs, by = "gene") 

dfhm$variable <- gsub("[0-9]", "", dfhm$variable)


# ggplot(dfhm[dfhm$variable %in% c("C", "G"),], aes(x = gene, y = value, fill = variable)) +
#   geom_boxplot() +
#   stat_compare_means(aes(group = variable),
#                      method = "anova",
#                      label =  "p.signif")+
#   facet_grid(~cat, scales = "free_x", space = "free_x")


# heat maps ----

deGs <- strsplit("Tbx21/ Ifngr1, Il18r1, Il18rap, Il12rb2/ Ifng", ", |/ ")[[1]]

deGs <- data.frame(gene = deGs, cat = c(rep("Transcription factor", 1), 
                                        rep("Receptors", 4),
                                        rep("Cytokiunes", 1)))

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- as.data.frame(dfhm) %>% 
  rownames_to_column("gene") %>% 
  reshape::melt(id = c("gene")) %>% 
  merge(deGs, by = "gene") 

dfhm$variable <- gsub("[0-9]", "", dfhm$variable)

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

p1 <- Heatmap(htmp,
              na_col = "black",
              rect_gp = gpar(col = "gray50", lwd = 0.8),
              row_dend_width = unit(2, "cm"),
              col = col_fun,
              row_names_gp = gpar(fontsize = 10),
              cluster_rows = T,
              show_row_dend = T,
              cluster_columns = F,
              column_split = c(1, rep(2, 3), 3, rep(4, 3)),
              row_split = c(rep(1, 1), rep(2, 4), rep(3, 1)),
              row_title = "Th1 Genes",
              column_names_side = "bottom",
              # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 10),
              column_title = NULL,
              show_column_names = F,
              name = "Normalized Gene Expression",
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#1B9E77","#D95F02", "#D95F02")),
                                                                  labels = c(
                                                                    expression(Uninfected-""), 
                                                                    expression(Infected-Foxp3^{Cre-YFP}),
                                                                    expression(Uninfected-""),
                                                                    expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                                  ),
                                                                  labels_gp = gpar(col = "white", fontsize = 9))),
              left_annotation = rowAnnotation(foo = c("Transcription factors", rep("Receptor", 4), "Cytokines"),
                                              col = list(foo = c("Transcription factors" = "gray30", "Receptor" = "brown", "Cytokines" = "blue")),
                                              show_annotation_name = c(foo = FALSE),
                                              annotation_legend_param = list(foo = list(title = "Gene Type", nrow = 1,
                                                                                        title_position = "topcenter",
                                                                                        direction = "horizontal"))),
              
              heatmap_legend_param = list(direction = "horizontal",
                                          title_position = "topcenter",
                                          legend_width = unit(7.5, "cm")))

deGs <- strsplit("Batf, Gata3/ Il4", ", |/ ")[[1]]
deGs <- data.frame(gene = deGs, cat = c(rep("Transcription factor", 2), 
                                        # rep("Receptors", 1),
                                        rep("Cytokiunes", 1)))

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- as.data.frame(dfhm) %>% 
  rownames_to_column("gene") %>% 
  reshape::melt(id = c("gene")) %>% 
  merge(deGs, by = "gene") 

dfhm$variable <- gsub("[0-9]", "", dfhm$variable)

dfhm <- AverageExpression(dat,features = deGs$gene, group.by = "reps", assays = "RNA")[[1]]

dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

p2 <- Heatmap(htmp,
              na_col = "black",
              rect_gp = gpar(col = "gray50", lwd = 0.8),
              row_dend_width = unit(2, "cm"),
              col = col_fun,
              row_names_gp = gpar(fontsize = 10),
              cluster_rows = T,
              show_row_dend = T,
              cluster_columns = F,
              column_split = c(1, rep(2, 3), 3, rep(4, 3)),
              row_split = c(rep(1, 2), rep(2, 1)),
              row_title = "Th2 Genes",
              column_names_side = "bottom",
              # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 10),
              column_title = NULL,
              show_column_names = F,
              name = "Normalized Gene Expression",
              # top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#1B9E77","#D95F02", "#D95F02")),
              #                                                     labels = c(
              #                                                       expression(Uninfected-""),
              #                                                       expression(Infected-Foxp3^{Cre-YFP}),
              #                                                       expression(Uninfected-""),
              #                                                       expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
              #                                                     ),
              #                                                     labels_gp = gpar(col = "white", fontsize = 9))),
              left_annotation = rowAnnotation(foo = c("Transcription factors", "Transcription factors", "Cytokines"),
                                              col = list(foo = c("Transcription factors" = "gray30", "Cytokines" = "blue")),
                                              show_annotation_name = c(foo = FALSE),
                                              annotation_legend_param = list(foo = list(title = "Gene Type", nrow = 1,
                                                                                        title_position = "topcenter",
                                                                                        direction = "horizontal"))),
              heatmap_legend_param = list(direction = "horizontal",
                                          title_position = "topcenter",
                                          legend_width = unit(15, "cm")))
png("3.Output/treg/selectedGenes_C1-3vsG1-Th1nTh2.png", height = 5, width = 9, units = "in", res = 300)
draw(p1 %v% p2, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()
# old ----
# rownames(htmp)[1] <- "Tbx21*"
# rownames(htmp)[5] <- "Ifngr1**"
# p1 <- draw(Heatmap(htmp,
#                    na_col = "black",
#                    rect_gp = gpar(col = "gray50", lwd = 0.8),
#                    row_dend_width = unit(2, "cm"),
#                    col = col_fun,
#                    row_names_gp = gpar(fontsize = 10),
#                    cluster_rows = F,
#                    show_row_dend = F,
#                    cluster_columns = F,
#                    column_split = c(1, rep(2, 3), 3, rep(4, 3)),
#                    row_split = c(rep(1, 1), rep(2, 4), rep(3, 1)),
#                    row_title = NULL,
#                    column_names_side = "bottom",
#                    # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
#                    column_names_rot = 0,
#                    column_names_gp = gpar(fontsize = 10),
#                    column_title = NULL,
#                    show_column_names = F,
#                    name = "Normalized Gene Expression",
#                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#1B9E77","#D95F02", "#D95F02")),
#                                                                        labels = c(
#                                                                          expression(Uninfected-""), 
#                                                                          expression(Infected-Foxp3^{Cre-YFP}),
#                                                                          expression(Uninfected-""),
#                                                                          expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
#                                                                        ),
#                                                                        labels_gp = gpar(col = "white", fontsize = 9))),
#                    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("gray30", "brown","blue")),
#                                                                     labels = c(
#                                                                       # "Transcription factors", 
#                                                                       # "Receptors",
#                                                                       # "Cytokines"
#                                                                     ),
#                                                                     labels_gp = gpar(col = "white", fontsize = 10))),
#                    heatmap_legend_param = list(direction = "horizontal",
#                                                title_position = "topcenter",
#                                                legend_width = unit(7.5, "cm"))), 
#            heatmap_legend_side = "bottom")
# 
# 
# # rownames(htmp)[1] <- "Batf*"
# # rownames(htmp)[6] <- "Il4*"
# p2 <- draw(Heatmap(htmp,
#              na_col = "black",
#              rect_gp = gpar(col = "gray50", lwd = 0.8),
#              row_dend_width = unit(2, "cm"),
#              col = col_fun,
#              row_names_gp = gpar(fontsize = 10),
#              cluster_rows = F,
#              show_row_dend = T,
#              cluster_columns = F,
#              column_split = c(1, rep(2, 3), 3, rep(4, 3)),
#              row_split = c(rep(1, 2), rep(2, 1)),
#              row_title = NULL,
#              column_names_side = "bottom",
#              # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
#              column_names_rot = 0,
#              column_names_gp = gpar(fontsize = 10),
#              column_title = NULL,
#              show_column_names = F,
#              name = "Normalized Gene Expression",
#              # top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#1B9E77","#D95F02", "#D95F02")),
#              #                                                     labels = c(
#              #                                                       expression(Uninfected-""), 
#              #                                                       expression(Infected-Foxp3^{Cre-YFP}),
#              #                                                       expression(Uninfected-""),
#              #                                                       expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
#              #                                                     ),
#              #                                                     labels_gp = gpar(col = "white", fontsize = 9))),
#              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("gray30", "brown")),
#                                                                  labels = c(
#                                                                    # "Transcription factors",
#                                                                    # "Cytokines", 
#                                                                    # "Receptors"
#                                                                  ),
#                                                                  labels_gp = gpar(col = "white", fontsize = 10))),
#              heatmap_legend_param = list(direction = "horizontal",
#                                          title_position = "topcenter",
#                                          legend_width = unit(7.5, "cm"))), heatmap_legend_side = "bottom")
# 
# png("3.Output/treg/DEgenes_C1-3vsG1-Th1.png", height = 6, width = 8, units = "in", res = 300)
# p1
# dev.off()
# 
# png("3.Output/treg/DEgenes_C1-3vsG1-Th2.png", height = 6, width = 8, units = "in", res = 300)
# p2
# dev.off()
# 
# p1 %v% p2
