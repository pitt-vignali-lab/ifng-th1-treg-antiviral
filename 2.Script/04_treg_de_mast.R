# ----------------------------------------------------------------------------
# Analysis | Differential expression between Ifngr1-deficient and control Tregs
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/3.1.deMast.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/deMAST_*.csv
#            3.Output/treg/deGsForVolcanoPlot.csv
# Upstream : 02_treg_cluster_labeling.R
# Note     : Seurat FindMarkers with the MAST test, logFC threshold 0.05
# ----------------------------------------------------------------------------

# 1. Data initialization and labelling ----
library(Seurat)
library(scales)
library(tidyverse)

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
               "Metabolically hyperactive", # 2
               "Proliferating", # 3
               "Gene transcription and protein stability" # 4
)

names(clusterId) <- levels(dat)

datLab <- RenameIdents(dat, clusterId)

DimPlot(datLab, label = T) + NoLegend()

# 2. DE analaysis ----

for (i in clusterId) {
  datSub <- subset(datLab, idents = i)
  Idents(datSub) <- "samples"
  
  FindMarkers(datSub, ident.1 = "G", ident.2 = "C",
              test.use = "MAST", logfc.threshold =  0.05) %>% 
    filter(p_val_adj < 0.05) %>% 
    rownames_to_column(var = "gene") %>% 
    write.csv(paste0("3.Output/treg/deMAST_", i, ".csv"))
}

Idents(datLab) <- "samples"

FindMarkers(datLab, ident.1 = "G", ident.2 = "C",
            test.use = "MAST",  logfc.threshold =  0.01)%>% 
  filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene")%>% 
  write.csv(paste0("3.Output/treg/deMAST_Overall.csv"))

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
           # write.csv(paste0("3.Output/treg/deMAST_", i, ".csv"))
         )
}

Idents(datLab) <- "samples"

deOverall <- FindMarkers(datLab, ident.1 = "G", ident.2 = "C",
                         test.use = "MAST",  
                         logfc.threshold =  0, 
                         min.pct = 0)%>% # was logfc = 0.05
  # filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene")# %>% 
  # write.csv(paste0("3.Output/treg/deMAST_Overall.csv"))
write.csv(deOverall, "3.Output/treg/deGsForVolcanoPlot.csv")

## volcano ----
library(ggrepel)

de <- read.csv("3.Output/treg/deGsForVolcanoPlot.csv")
volDf <- data.frame()
for (j in ls()[grep("de", ls())][1]) {
  tDf <- eval(parse(text = paste0("`", j, "`")))
  tDf$cellType <- gsub("de", "", j)
  volDf <- rbind(volDf, tDf)
}


volDf$cellType <- factor(volDf$cellType, level = c("Overall", "Central", "Effector",
                                                       "Gene transcription and protein stability",
                                                       "Metabolically hyperactive",
                                                       "Proliferating"))
gsHl <- c("Stat1", "Gata3", "Il7r", "Ifngr1", "Cxcr3", "Klrg1", 
          "Tbx21", "AW112010", "Ly6c2")

gsAll <- unique(c(gsHl,
                  head(volDf[order(volDf$p_val) & volDf$avg_log2FC > 0, ]$gene),
                  head(volDf[order(volDf$p_val) & volDf$avg_log2FC < 0, ]$gene)))



png("3.Output/treg/deVolcanoOverall.png", width = 8.5, height = 8, units = "in", res = 300)
ggplot(volDf, 
       aes(x = avg_log2FC, y = -log(p_val), 
           col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                        ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN", NA)))) + 
  geom_point(size = 1.8, 
             alpha = ifelse(volDf$gene %in% gsHl, 0.8, 0.25), 
             # shape = ifelse(volDf$gene %in% gsHl, 19, 1)
             ) +
  xlim(-1, 1) +
  geom_text_repel(data = volDf[volDf$gene %in% gsAll,], 
    aes(x = avg_log2FC, y = -log(p_val), label = gene,
        col = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DN", NA))),
    fontface = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, "bold", "plain"),
    alpha = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 1, 0.7),
    size = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 5, 4.5),
    show.legend  = FALSE, 
    nudge_x = 0.05) +
  scale_color_manual(name = "Status",
                     breaks = c("DN", "UP"),
                     values=c("#FF0000", "#7030A0"),
                     label = c(expression(paste(UP," ","in"," ",Foxp3^{Cre-YFP})),
                               expression(paste(UP," ","in"," ",Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
  # facet_wrap(~cellType, ncol = 3) + 
  ylab("Significance (-logPval)") +
  xlab("Log fold change") +
  # xlim(-2, 2) +
  theme_bw() +
  # guides(col=guide_legend(nrow=2,byrow=TRUE),
  #        text = NULL) +
  theme(text = element_text(size=15),
        legend.position = "bottom",
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
png("3.Output/treg/Venn6cat.png", width = 6, height = 6, units = "in", res = 300)
venn(venLs[c(1, 3, 6, 2, 4, 7)], zcolor = "style")
dev.off()

png("3.Output/treg/Venn9cat.png", width = 6, height = 6, units = "in", res = 300)
venn(venLs[c(1, 3, 6, 8, 5, 4, 2, 7, 9)], zcolor = "style")
dev.off()

## heatmap ----
library(ComplexHeatmap)
library(BBmisc)
library(scales)   
### Overall ----

dfhm <- AverageExpression(datLab,features = deOverall[deOverall$p_val_adj<0.05,]$gene, 
                          group.by = "reps", assays = "RNA")[[1]]
dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

png("3.Output/treg/deHeatmapOverall.png", width = 9, height = 16, units = "in", res = 300)
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
        column_split = c(1, rep(2, 3), 3, rep(4, 3)),
        # row_split = c(rep(1, 1), rep(2, 4), rep(3, 1)),
        # row_title = "Th1 Genes",
        column_names_side = "bottom",
        # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
        column_names_rot = 0,
        column_names_gp = gpar(fontsize = 10),
        column_title = NULL,
        show_column_names = F,
        name = "Normalized Gene Expression",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#ff0000", "#ff0000","#7030a0", "#7030a0")),
                                                            labels = c(
                                                              expression(Uninfected-""), 
                                                              expression(Infected-Foxp3^{Cre-YFP}),
                                                              expression(Uninfected-""),
                                                              expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                            ),
                                                            labels_gp = gpar(col = "white", fontsize = 9))),
        heatmap_legend_param = list(direction = "horizontal",
                                    title_position = "topcenter",
                                    legend_width = unit(7.5, "cm"))),
heatmap_legend_side = "bottom")
dev.off()

### seperate ----
dfhm <- data.frame()
gs <- c()
for(j in clusterId[-5]){
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
dfhm <- dfhm[,c(7, 1, 2, 3, 8, 4, 5, 6)]

htmp <- normalize(as.matrix(dfhm), range = c(-2, 2))

colnames(htmp) <- colnames(dfhm)

png("3.Output/treg/deHeatmapSep.png", width = 9, height = 16, units = "in", res = 300)
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
        column_split = c(1, rep(2, 3), 3, rep(4, 3)),
        row_split = rep(1:4, gs),
        # column_names_side = "bottom",
        show_row_name = F,
        row_title = NULL,
        # row_title = gsub("DE_", "", ls()[grep("DE", ls())][i]),
        # column_names_rot = 0,
        # column_names_gp = gpar(fontsize = 10),
        column_title = NULL,
        show_column_name = F,
        name = "Normalized Gene Expression",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#1B9E77", "#1B9E77","#D95F02", "#D95F02")),
                                                            labels = c(
                                                              expression(Uninfected-""), 
                                                              expression(Infected-Foxp3^{Cre-YFP}),
                                                              expression(Uninfected-""),
                                                              expression(paste(Infected-Ifngr1^{L/L},Foxp3^{Cre-YFP}))
                                                            ),
                                                            labels_gp = gpar(col = "white", fontsize = 9))),
        left_annotation = rowAnnotation(foo = rep(clusterId[-5], gs),
                                        col = list(foo = c("Effector" = hue_pal()(5)[1], 
                                                           "Central" = hue_pal()(5)[2],
                                                           "Metabolically hyperactive" = hue_pal()(5)[3],
                                                           "Proliferating" = hue_pal()(5)[4])),
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
