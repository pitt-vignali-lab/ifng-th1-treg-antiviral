# ----------------------------------------------------------------------------
# Analysis | Differential expression in tetramer-positive CD8 T cells
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/3.1.deAnalysis.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/deMAST_Overall.csv
#            3.Output/cd8/deVolcano.png
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
require(scales)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")
DefaultAssay(dat) <- "RNA"
Idents(dat) <- "RNA_snn_res.0.3"

datSub <- subset(dat, idents = 8, invert = TRUE)

clusterId <- c("Activated", # 0
               "Cytotoxic", # 1
               "Progenitor exhausted", # 2
               "Proliferating", # 3
               "Exhausted", # 4
               "Memory", # 5
               "IFN signaling", # 6
               "Metabolically hyperactive" # 7
)

names(clusterId) <- levels(datSub)

datLab <- RenameIdents(datSub, clusterId)

# 2 DE ----
for (i in clusterId) {
  datSub <- subset(datLab, idents = i)
  Idents(datSub) <- "samples"
  
  FindMarkers(datSub, ident.1 = "G", ident.2 = "C",
              test.use = "MAST", logfc.threshold =  0.05) %>% 
    filter(p_val_adj < 0.05) %>% 
    rownames_to_column(var = "gene") %>% 
    write.csv(paste0("3.Output/cd8/deMAST_", i, ".csv"))
}

Idents(datLab) <- "samples"

FindMarkers(datLab, ident.1 = "G", ident.2 = "C",
            test.use = "MAST",  logfc.threshold =  0.01)%>% 
  filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene")%>% 
  write.csv(paste0("3.Output/cd8/deMAST_Overall.csv"))

# 3. visualization ----
datLab <- RenameIdents(dat, clusterId)
for (i in clusterId) {
  datSub <- subset(datLab, idents = i)
  Idents(datSub) <- "samples"
  assign(paste0("de", i), 
         FindMarkers(datSub, ident.1 = "G", ident.2 = "C",
                     test.use = "MAST", logfc.threshold =  0.05) %>% 
           #filter(p_val_adj < 0.05) %>% 
           rownames_to_column(var = "gene") #%>% 
         # write.csv(paste0("3.Output/cd8/deMAST_", i, ".csv"))
  )
}

Idents(datLab) <- "samples"

deOverall <- FindMarkers(datLab, ident.1 = "G", ident.2 = "C",
                         test.use = "MAST",  logfc.threshold =  0.01)%>% 
  # filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene")# %>% 
# write.csv(paste0("3.Output/cd8/deMAST_Overall.csv"))

## volcano ----
library(ggrepel)

volDf <- data.frame()
for (j in ls()[grep("de", ls())]) {
  tDf <- eval(parse(text = paste0("`", j, "`")))
  tDf$cellType <- gsub("de", "", j)
  volDf <- rbind(volDf, tDf)
}

volDf$cellType <- factor(volDf$cellType, level = c("Overall", 
                                                   "Progenitor exhausted", # 2
                                                   "Activated", # 0
                                                   "Exhausted", # 4
                                                   "Cytotoxic", # 1
                                                   "Proliferating", # 3
                                                   "IFN signaling", # 6
                                                   "Memory", # 5
                                                   "Metabolically hyperactive"))

png("3.Output/cd8/deVolcano.png", width = 10, height = 10, units = "in", res = 300)
ggplot(volDf, 
       aes(x = avg_log2FC, y = -log(p_val), 
           col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                        ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN", NA)))) + 
  geom_point(size = 0.5) +
  geom_text_repel(data = volDf[volDf$p_val_adj < 0.05,], 
                  aes(x = avg_log2FC, y = -log(p_val), label = gene,
                      col = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DN", NA))),
                  size = 2.5) +
  scale_color_manual(name = "Status",
                     breaks = c("UP", "DN"),
                     values=c("#D95F02", "#1B9E77"),
                     label = c("significantly up-regulated genes in the mutant group", 
                               "significantly up-regulated genes in the control group")) +
  facet_wrap(~cellType, ncol = 3) + 
  ylab("Significance (-logPval)") +
  xlab("Log fold change") +
  xlim(-2, 2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.key=element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"))
dev.off()

## venn ----
library(VennDiagram)

venLs <- list()
for (j in ls()[grep("de", ls())]) {
  tLs <- eval(parse(text = paste0("`", j, "`")))
  tLs <- list(tLs[tLs$p_val_adj < 0.05 & tLs$avg_log2FC > 0,]$gene,
              tLs[tLs$p_val_adj < 0.05 & tLs$avg_log2FC < 0,]$gene)
  names(tLs) <- paste0(gsub("de", "", j), c("_UP", "_DN"))
  venLs <- c(venLs, tLs)
}

venLs <- venLs[lengths(venLs) != 0]

VennDiagram <- venn.diagram(x = venLs[1:4], filename = NULL)
cowplot::plot_grid(VennDiagram)

genes <- paste0("gene",1:1000)
gene_list <- list(A = sample(genes,100),
                  B = sample(genes,200),
                  C = sample(genes,300),
                  D = sample(genes,200))


library(venn)
png("3.Output/cd8/Venn6cat.png", width = 6, height = 6, units = "in", res = 300)
venn(venLs[c(1, 3, 6, 2, 4, 7)], zcolor = "style")
dev.off()

png("3.Output/cd8/Venn9cat.png", width = 6, height = 6, units = "in", res = 300)
venn(venLs[c(1, 3, 6, 8, 5, 4, 2, 7, 9)], zcolor = "style")
dev.off()

## heatmap ----
library(ComplexHeatmap)
library(BBmisc)
library(scales)   
### Overall ----

dfhm <- AverageExpression(datLab,features = deOverall[deOverall$p_val_adj<0.05,]$gene, 
                          group.by = "reps", assays = "RNA")[[1]]
# dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

png("3.Output/cd8/deHeatmapOverall.png", width = 9, height = 16, units = "in", res = 300)
draw(
  Heatmap(htmp,
          na_col = "black",
          rect_gp = gpar(col = "gray50", lwd = 0.8),
          row_dend_width = unit(2, "cm"),
          col = col_fun,
          row_names_gp = gpar(fontsize = 10),
          cluster_rows = T,
          show_row_dend = T,
          cluster_columns = F,
          column_split = c(rep(1, 3), rep(2, 3)),
          # row_split = c(rep(1, 1), rep(2, 4), rep(3, 1)),
          # row_title = "Th1 Genes",
          column_names_side = "bottom",
          # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
          column_names_rot = 0,
          column_names_gp = gpar(fontsize = 10),
          column_title = NULL,
          show_column_names = F,
          name = "Normalized Gene Expression",
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#D95F02")),
                                                              labels = c(
                                                                # expression(Uninfected-""), 
                                                                expression(Infected-Foxp3^{Cre-YFP}),
                                                                # expression(Uninfected-""),
                                                                expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                              ),
                                                              labels_gp = gpar(col = "white", fontsize = 9))),
          heatmap_legend_param = list(direction = "horizontal",
                                      title_position = "topcenter",
                                      legend_width = unit(7.5, "cm"))),
  heatmap_legend_side = "bottom")
dev.off()

### seperate ----
datLab <- RenameIdents(dat, clusterId)
dfhm <- data.frame()
gs <- c()
for(j in clusterId[-c(6:8)]){
  datSub <- subset(datLab, idents = j)
  des <- eval(parse(text = paste0("`de", j, "`")))
  tdf <- AverageExpression(datSub,features = des[des$p_val_adj<0.05,]$gene, 
                           group.by = "reps", assays = "RNA")[[1]] 
  rownames(tdf) <- des[des$p_val_adj<0.05,]$gene
  tdf<- as.data.frame(tdf) %>% 
    rownames_to_column("gene") 
  gs <- c(gs, nrow(tdf))
  dfhm <- rbind(dfhm, tdf)
}

labels <- dfhm$gene
dfhm <- dfhm[-1]
# dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

png("3.Output/cd8/deHeatmapSep.png", width = 9, height = 12, units = "in", res = 300)
draw(
  Heatmap(htmp,
          na_col = "black",
          rect_gp = gpar(col = "gray50", lwd = 0.8),
          row_dend_width = unit(2, "cm"),
          col = col_fun,
          # row_names_gp = gpar(fontsize = 0),
          cluster_rows = T,
          show_row_dend = T,
          cluster_columns = F,
          column_split = c(rep(1, 3), rep(2, 3)),
          row_split = rep(1:5, gs),
          # column_names_side = "bottom",
          show_row_name = F,
          row_title = NULL,
          # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
          # column_names_rot = 0,
          # column_names_gp = gpar(fontsize = 10),
          column_title = NULL,
          show_column_name = F,
          name = "Normalized Gene Expression",
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#D95F02")),
                                                              labels = c(
                                                                # expression(Uninfected-""), 
                                                                expression(Infected-Foxp3^{Cre-YFP}),
                                                                # expression(Uninfected-""),
                                                                expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                              ),
                                                              labels_gp = gpar(col = "white", fontsize = 9))),
          left_annotation = rowAnnotation(foo = rep(clusterId[-c(6:8)], gs),
                                          col = list(foo = c("Progenitor exhausted"= hue_pal()(8)[1], # 2
                                                             "Activated"= hue_pal()(8)[2], # 0
                                                             "Exhausted"= hue_pal()(8)[3], # 4
                                                             "Cytotoxic"= hue_pal()(8)[4], # 1
                                                             "Proliferating"= hue_pal()(8)[5])),
                                          show_annotation_name = c(foo = FALSE),
                                          annotation_legend_param = list(foo = list(title = "Cell type", nrow = 1,
                                                                                    title_position = "topcenter",
                                                                                    direction = "horizontal"))),
          heatmap_legend_param = list(direction = "horizontal",
                                      title_position = "topcenter",
                                      legend_width = unit(7.5, "cm")))+ 
    rowAnnotation(labels = anno_text(labels, which = "row"), 
                  width = unit(3, "cm")),
  heatmap_legend_side = "bottom", 
  annotation_legend_side = "bottom")
dev.off()

