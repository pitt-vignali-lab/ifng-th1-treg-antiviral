# ----------------------------------------------------------------------------
# Revision | Pseudobulk heat maps and volcano plots
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/5.2.heatmaps&volcanaol.R
# Input    : 3.Output/revision/4.PseudoBulk/
# Output   : 3.Output/revision/4.PseudoBulk/volcanol_*.png
# Upstream : 35_rev_pseudobulk.R
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/5.PseudoBulk")){
  dir.create("3.Output/revision/5.PseudoBulk")
}
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(extrafont)
library(GetoptLong)  # for the function qq()

# data read in ----

dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")

ctMat <- AverageExpression(dat, 
                           group.by = "reps", 
                           alot  = "data",
                           assays = "RNA")[[1]]

ctMat1 <- normalize(as.matrix(ctMat[,c(1:3, 6:8, 11, 14)]), 
                    range = c(-2, 2))
colnames(ctMat1) <- colnames(ctMat[,c(1:3, 6:8, 11, 14)])

ctMat2 <- normalize(as.matrix(ctMat[,c(4, 5, 9, 10, 12, 13, 15, 16)]), 
                    range = c(-2, 2))
colnames(ctMat2) <- colnames(ctMat[,c(4, 5, 9, 10, 12, 13, 15, 16)])

ctMat <- cbind(ctMat1, ctMat2)

ctMat <- ctMat[,c(paste0("UC", 1:3),
                  paste0("C", 1:5),
                  paste0("UG", 1:3),
                  paste0("G", 1:5))]


deG <- read.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv")

deUG <- read.csv("3.Output/revision/5.PseudoBulk/deDfUGvsUC.csv")


# heatmap overall settings ----

group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

split = c(rep(1, 3), rep(2, 5), rep(3, 3), rep(4, 5))

ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = c("#fb0707", "#fb0707",
                                      "#7434a4", "#7434a4")),
                   labels = c("Uninfected", "Infected", 
                              "Uninfected", "Infected" ),
                   labels_gp = gpar(col = "white", 
                                    fontsize = 12,
                                    fontface = 2)))
# heatmap 1 (seletcted genes) ----
htGenes <- stringr::str_to_title(strsplit("Ifngr1, Ifngr2, Il12rb2, Il12rb1, TIGIT, Ikzf2, ICOS, STAT5a, Ikzf4, FOXO1, Foxo3, CTLA4, TNFRSF18, Nrp1, ENTPD1, NT5E, LAG3, PDCD1, CD274, IL10, Ebi3, Il12b, TGFB1, IL2RA, CXCR3, Il16, Ccl5, Ly6a", ", ")[[1]])

htGenes <- stringr::str_to_title(strsplit("Il16, Ccr5, Casp3, Cd4, Ccl5", ", ")[[1]])

htGenes <- deG %>% 
  filter(feature %in% htGenes) %>%
  pull(feature)

htmp <- normalize(ctMat[htGenes, ], range = c(-2, 2))


rownames(htmp) <- paste0(rownames(htmp), 
                         " (padj: ", 
                         format(sapply(rownames(htmp), 
                                       function(x) deG[deG$feature == x,]$padj), 
                                digits = 3, scientific = T), 
                         ")")

fontCol <- ifelse(sapply(htGenes, 
                         function(x) {
                           deG[deG$feature == x,]$logFC}) == 0, 
                  "gray50", 
                  ifelse(sapply(htGenes, 
                         function(x) {
                           deG[deG$feature == x,]$logFC}) > 0, 
                  "#7434a4", "#fb0707"))


png("3.Output/revision/5.PseudoBulk/selectedGCytokineHeatmap.png",
    width = 8.5, height = 3, units = "in", res = 300)
draw(Heatmap(
  htmp,
  na_col = "black",
  rect_gp = gpar(col = "gray50", lwd = 0.8),
  row_dend_width = unit(2, "cm"),
  col = col_fun, 
  row_names_gp = gpar(fontsize = 12, 
                      fontface = 2,
                      col = fontCol),
  cluster_rows = T,
  show_row_dend = T,
  cluster_columns = F,
  
  column_title = NULL,
  row_names_side = "left",
  show_column_names = F,
  name = "Normalized Gene Expression",
  column_split = split, 
  top_annotation = ha,
  heatmap_legend_param = list(direction = "horizontal",
                              title_position = "topcenter",
                              legend_width = unit(7.5, "cm"))),
  heatmap_legend_side = "bottom")

group_block_anno(1:2, "empty", gp = gpar(fill = "#fb0707"), 
                 label = bquote(bold(Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))

group_block_anno(3:4, "empty", gp = gpar(fill = "#7434a4"), 
                 label = bquote(bold(Ifngr1^{L/L}~Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))


dev.off()

# heatmao 2 ???Foxp3 induced and repressed genes)----
## targeted genes ----
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

## plotting ----
htGenes <- data.frame(gene = c(gIndu, gRepr),
                   cat = c(rep("Foxp3 Induced", length(gIndu)),
                           rep("Foxp3 Repressed", length(gRepr))))

dfhm <- ctMat[htGenes$gene,] 

# dfhm <- dfhm[,c(11:13, 1:5, 14:16, 6:10)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

rownames(htmp) <- paste0(rownames(htmp), 
                         " (padj: ", 
                         format(sapply(rownames(htmp), 
                                       function(x) deG[deG$feature == x,]$padj), 
                                digits = 3, scientific = T), 
                         ")")

fontCol <- ifelse(sapply(htGenes$gene, 
                         function(x) {
                           deG[deG$feature == x,]$logFC}) == 0, 
                  "gray50", 
                  ifelse(sapply(htGenes$gene, 
                                function(x) {
                                  deG[deG$feature == x,]$logFC}) > 0, 
                         "#7434a4", "#fb0707"))

sidAnno <- rowAnnotation(foo = anno_block(gp = gpar(fill = c("#7C4935", "#457373")),
                                          labels = c("Foxp3 induced genes", "Foxp3 repressed genes"),
                                          labels_gp = gpar(col = "white", 
                                                           fontsize = 10,
                                                           fontface = 2)))

png("3.Output/revision/5.PseudoBulk/foxp3Genesheat.png", height = 16, width = 9, units = "in", res = 300)

draw(Heatmap(
  htmp,
  na_col = "black",
  rect_gp = gpar(col = "gray50", lwd = 0.8),
  row_dend_width = unit(2, "cm"),
  col = col_fun, 
  row_names_gp = gpar(fontsize = 12, 
                      fontface = 2,
                      col = fontCol),
  
  left_annotation = sidAnno,
  row_split = c(rep("Indu", length(gIndu)),
                rep("Repres", length(gRepr))),
  cluster_rows = T,
  show_row_dend = T,
  cluster_columns = F,
  column_title = NULL,
  row_title = NULL,
  row_names_side = "left",
  show_column_names = F,
  name = "Normalized Gene Expression",
  column_split = split, 
  top_annotation = ha,
  heatmap_legend_param = list(direction = "horizontal",
                              title_position = "topcenter",
                              legend_width = unit(7.5, "cm"))),
  heatmap_legend_side = "bottom")

group_block_anno(1:2, "empty", gp = gpar(fill = "#fb0707"), 
                 label = bquote(bold(Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))

group_block_anno(3:4, "empty", gp = gpar(fill = "#7434a4"), 
                 label = bquote(bold(Ifngr1^{L/L}~Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))
dev.off()

# heatmap3 (Th1 and Th2 genes)----

htGenes <- c(strsplit("Tbx21/ Ifngr1, Il18r1, Il18rap, Il12rb2/ Ifng", ", |/ ")[[1]],
             strsplit("Batf, Gata3/ Il4", ", |/ ")[[1]])

dfhm <- ctMat[htGenes,] 

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

rownames(htmp) <- paste0(rownames(htmp), 
                         " (padj: ", 
                         format(sapply(rownames(htmp), 
                                       function(x) deG[deG$feature == x,]$padj), 
                                digits = 3, scientific = T), 
                         ")")

fontCol <- ifelse(sapply(htGenes, 
                         function(x) {
                           deG[deG$feature == x,]$logFC}) == 0, 
                  "gray50", 
                  ifelse(sapply(htGenes, 
                                function(x) {
                                  deG[deG$feature == x,]$logFC}) > 0, 
                         "#7434a4", "#fb0707"))


png("3.Output/revision/5.PseudoBulk/th1th2genes.png",
    width = 9, height = 4.5, units = "in", res = 300)
draw(Heatmap(htmp,
             na_col = "black",
             rect_gp = gpar(col = "gray50", lwd = 0.8),
             row_dend_width = unit(2, "cm"),
             col = col_fun,
             row_names_gp = gpar(fontsize = 12, 
                                 fontface = 2,
                                 col = fontCol),
             cluster_rows = T,
             show_row_dend = T,
             cluster_columns = F,
             column_split = split,
             row_split = data.frame(s1 = c(rep("Th1 Genes", 6), rep("Th2 Genes", 3)),
                                    s2 = c(rep(1, 1), rep(2, 4), rep(3, 1), rep(4, 2), rep(5, 1))),
             row_title = c("", "Th1 Genes", "", "", "Th2 Genes"),
             # cluster_row_slices = FALSE,
             column_names_side = "bottom",
             column_names_rot = 0,
             column_names_gp = gpar(fontsize = 10),
             column_title = NULL,
             show_column_names = F,
             name = "Normalized Gene Expression",
             top_annotation = ha,
             left_annotation = rowAnnotation(foo = c("Transcription factors", 
                                                     rep("Receptor", 4), 
                                                     "Cytokines",
                                                     rep("Transcription factors", 2), 
                                                     "Cytokines"),
                                             col = list(foo = c("Transcription factors" = "gray30", 
                                                                "Receptor" = "brown", 
                                                                "Cytokines" = "blue")),
                                             show_annotation_name = c(foo = FALSE),
                                             annotation_legend_param = list(foo = list(title = "Gene Type", nrow = 1,
                                                                                       title_position = "topcenter",
                                                                                       direction = "horizontal"))),
             
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = unit(7.5, "cm"))), 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")

group_block_anno(1:2, "empty", gp = gpar(fill = "#fb0707"), 
                 label = bquote(bold(Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))

group_block_anno(3:4, "empty", gp = gpar(fill = "#7434a4"), 
                 label = bquote(bold(Ifngr1^{L/L}~Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))
dev.off()

# heatmap4 (top ranked genes) ----
htGenes <- deG |> 
  arrange(padj) |> 
  head(50) |> 
  pull(feature)

dfhm <- ctMat[htGenes,]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2),
                  method = "standardize")

rownames(htmp) <- paste0(rownames(htmp), 
                         " (padj: ", 
                         format(sapply(rownames(htmp), 
                                       function(x) deG[deG$feature == x,]$padj), 
                                digits = 3, scientific = T), 
                         ")")

fontCol <- ifelse(is.na(sapply(htGenes, 
                               function(x) {
                                 deG[deG$feature == x,]$logFC})), 
                  "gray50", 
                  ifelse(sapply(htGenes, 
                                function(x) {
                                  deG[deG$feature == x,]$logFC}) > 0, 
                         "#7434a4", "#fb0707"))

png("3.Output/revision/5.PseudoBulk/top50RankedGenes.png",
    width = 9, height = 12, units = "in", res = 300)
draw(Heatmap(
  htmp,
  na_col = "black",
  rect_gp = gpar(col = "gray50", lwd = 0.8),
  row_dend_width = unit(2, "cm"),
  col = col_fun, 
  row_names_gp = gpar(fontsize = 12, 
                      fontface = 2,
                      col = fontCol),
  cluster_rows = T,
  show_row_dend = T,
  cluster_columns = F,
  
  column_title = NULL,
  row_names_side = "left",
  show_column_names = F,
  name = "Normalized Gene Expression",
  column_split = split, 
  top_annotation = ha,
  heatmap_legend_param = list(direction = "horizontal",
                              title_position = "topcenter",
                              legend_width = unit(7.5, "cm"))),
  heatmap_legend_side = "bottom")

group_block_anno(1:2, "empty", gp = gpar(fill = "#fb0707"), 
                 label = bquote(bold(Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))

group_block_anno(3:4, "empty", gp = gpar(fill = "#7434a4"), 
                 label = bquote(bold(Ifngr1^{L/L}~Foxp3^{Cre-YFP})),
                 label_gp = gpar(col = "white", 
                                 fontsize = 12))


dev.off()

# volcano plots ----
library(ggrepel)


# volDf <- read.csv("3.Output/treg/deGsForVolcanoPlot.csv")

gsHl <- c("Stat1", "Gata3", "Il7r", "Ifngr1", "Cxcr3", "Klrg1", 
          "Tbx21", "AW112010", "Ly6c2", "Nkg7", "Ifit3")

# gsAll <- unique(c(gsHl,
#                   head(volDf[order(volDf$padj) & volDf$logFC > 0, ]$gene),
#                   head(volDf[order(volDf$padj) & volDf$logFC < 0, ]$gene)))
# 
# 
# 
# png("3.Output/treg/2C.png", width = 7, height = 7, units = "in", res = 300)
# ggplot(volDf, 
#        aes(x = avg_log2FC, y = -log(p_val), 
#            col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
#                         ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN", NA)))) + 
#   geom_point(size = ifelse(volDf$gene %in% gsHl, 3, 2), 
#              alpha = ifelse(volDf$gene %in% gsHl, 0.8, 0.3),
#              # stroke = ifelse(volDf$gene %in% gsHl, 0, 1.5),
#              # shape = ifelse(volDf$gene %in% gsHl, 19, 1)
#   ) +
#   xlim(-1, 1) +
#   geom_text_repel(data = volDf[volDf$gene %in% gsAll,], 
#                   aes(x = avg_log2FC, y = -log(p_val), label = gene,
#                       col = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DN", NA))),
#                   fontface = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, "bold", "plain"),
#                   alpha = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 1, 0.7),
#                   size = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 6, 5),
#                   show.legend  = FALSE) +
#   scale_color_manual(name = "Status",
#                      breaks = c("DN", "UP"),
#                      values=c("#FF0000", "#7434a4"),
#                      label = c(expression(paste(UP," ","in"," ",Foxp3^{Cre-YFP})),
#                                expression(paste(UP," ","in"," ",Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
#   # facet_wrap(~cellType, ncol = 3) + 
#   ylab("Significance (-logPval)") +
#   xlab("Log fold change") +
#   # xlim(-2, 2) +
#   theme_bw() +
#   # guides(col=guide_legend(nrow=2,byrow=TRUE),
#   #        text = NULL) +
#   theme(text = element_text(family = "Arial", size = 18, face = "bold"),
#         panel.border = element_rect(size = 2),
#         panel.grid = element_blank(),
#         legend.position = "bottom",
#         legend.box = "vertical",
#         legend.key=element_blank(),
#         strip.background = element_rect(colour = "white", fill = "white"))
# dev.off()

gsHl <- unique(c(deG |> 
  # arrange(padj) |> 
  mutate(grp = logFC > 0) |> 
  group_by(grp) |> 
  top_n(n = 10, wt = -padj) |> 
  pull(feature), 
  c("Stat1", "Gata3", "Il7r", "Ifngr1", "Cxcr3", "Klrg1", 
    "Tbx21", "AW112010", "Ly6c2", "Nkg7", "Ifit3")))


ggplot(deG, 
       aes(x = logFC, y = -log(padj), 
           col = ifelse(padj < 0.05 & logFC > 0, "UP",
                        ifelse(padj < 0.05 & logFC < 0, "DN", NA)))) + 
  geom_point(size = ifelse(deG$feature %in% gsHl, 3, 2), 
             alpha = ifelse(deG$feature %in% gsHl, 0.8, 0.3),
             # stroke = ifelse(volDf$gene %in% gsHl, 0, 1.5),
             # shape = ifelse(volDf$gene %in% gsHl, 19, 1)
  ) +
  xlim(-max(abs(deG$logFC)), max(abs(deG$logFC))) +
  geom_text_repel(data = deG[deG$feature %in% gsHl,], 
                  aes(x = logFC, y = -log(padj), label = feature,
                      col = ifelse(logFC > 0, "UP", ifelse(logFC < 0, "DN", NA))),
                  max.overlaps = Inf,
                  fontface =  "bold",
                  # alpha = 1,
                  size = 6,
                  show.legend  = FALSE) +
  scale_color_manual(name = "Status",
                     breaks = c("DN", "UP"),
                     values=c("#FF0000", "#7434a4"),
                     label = c(expression(paste(UP," ","in"," ",infected~Foxp3^{Cre-YFP})),
                               expression(paste(UP," ","in"," ",infected~Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
  # facet_wrap(~cellType, ncol = 3) + 
  ylab("Significance (-logPadj)") +
  xlab("Log fold change") +
  # xlim(-2, 2) +
  theme_bw() +
  # guides(col=guide_legend(nrow=2,byrow=TRUE),
  #        text = NULL) +
  theme(text = element_text(family = "sans", size = 18, face = "bold"),
        panel.border = element_rect(size = 2),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key=element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"))

ggsave("3.Output/revision/5.PseudoBulk/volcanol_GvsC.pdf", 
       height = 7, width = 9, units = "in")



gsHl <- unique(c(deUG |> 
  # arrange(padj) |> 
  mutate(grp = logFC > 0) |> 
  group_by(grp) |> 
  top_n(n = 10, wt = -padj) |> 
  pull(feature) , 
c("Stat1", "Gata3", "Il7r", "Ifngr1", "Cxcr3", "Klrg1", 
  "Tbx21", "AW112010", "Ly6c2", "Nkg7", "Ifit3")))


ggplot(deUG, 
       aes(x = logFC, y = -log(padj), 
           col = ifelse(padj < 0.05 & logFC > 0, "UP",
                        ifelse(padj < 0.05 & logFC < 0, "DN", NA)))) + 
  geom_point(size = ifelse(deUG$feature %in% gsHl, 3, 2), 
             alpha = ifelse(deUG$feature %in% gsHl, 0.8, 0.3),
             # stroke = ifelse(volDf$gene %in% gsHl, 0, 1.5),
             # shape = ifelse(volDf$gene %in% gsHl, 19, 1)
  ) +
  geom_text_repel(data = deUG[deUG$feature %in% gsHl & deUG$padj < 0.05,], 
                  aes(x = logFC, y = -log(padj), label = feature,
                      col = ifelse(logFC > 0, "UP", ifelse(logFC < 0, "DN", NA))),
                  max.overlaps = Inf,
                  fontface =  "bold",
                  # alpha = 1,
                  size = 6,
                  show.legend  = FALSE) +
  scale_color_manual(name = "Status",
                     breaks = c("DN", "UP"),
                     values=c("#FF0000", "#7434a4"),
                     label = c(expression(paste(UP," ","in"," ",uninfected~Foxp3^{Cre-YFP})),
                               expression(paste(UP," ","in"," ",uninfected~Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
  # facet_wrap(~cellType, ncol = 3) + 
  ylab("Significance (-logPadj)") +
  xlab("Log fold change") +
  xlim(-max(abs(deUG$logFC)), max(abs(deUG$logFC))) +
  theme_bw() +
  # guides(col=guide_legend(nrow=2,byrow=TRUE),
  #        text = NULL) +
  theme(text = element_text(family = "sans", size = 18, face = "bold"),
        panel.border = element_rect(size = 2),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key=element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"))

ggsave("3.Output/revision/5.PseudoBulk/volcanol_UGvsUC.pdf", 
       height = 7, width = 9, units = "in")
