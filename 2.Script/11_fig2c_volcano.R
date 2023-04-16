# ----------------------------------------------------------------------------
# Fig 2c | Volcano plot of differentially expressed genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2C.R
# Input    : 3.Output/treg/deGsForVolcanoPlot.csv
# Output   : 3.Output/treg/2C.png
# Upstream : 04_treg_de_mast.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(extrafont)
library(ggrepel)

volDf <- read.csv("3.Output/treg/deGsForVolcanoPlot.csv")

gsHl <- c("Stat1", "Gata3", "Il7r", "Ifngr1", "Cxcr3", "Klrg1", 
          "Tbx21", "AW112010", "Ly6c2", "Nkg7", "Ifit3")

gsAll <- unique(c(gsHl,
                  head(volDf[order(volDf$p_val) & volDf$avg_log2FC > 0, ]$gene),
                  head(volDf[order(volDf$p_val) & volDf$avg_log2FC < 0, ]$gene)))



png("3.Output/treg/2C.png", width = 7, height = 7, units = "in", res = 300)
ggplot(volDf, 
       aes(x = avg_log2FC, y = -log(p_val), 
           col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                        ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN", NA)))) + 
  geom_point(size = ifelse(volDf$gene %in% gsHl, 3, 2), 
             alpha = ifelse(volDf$gene %in% gsHl, 0.8, 0.3),
             # stroke = ifelse(volDf$gene %in% gsHl, 0, 1.5),
             # shape = ifelse(volDf$gene %in% gsHl, 19, 1)
  ) +
  xlim(-1, 1) +
  geom_text_repel(data = volDf[volDf$gene %in% gsAll,], 
                  aes(x = avg_log2FC, y = -log(p_val), label = gene,
                      col = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DN", NA))),
                  fontface = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, "bold", "plain"),
                  alpha = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 1, 0.7),
                  size = ifelse(volDf[volDf$gene %in% gsAll,]$gene %in% gsHl, 6, 5),
                  show.legend  = FALSE) +
  scale_color_manual(name = "Status",
                     breaks = c("DN", "UP"),
                     values=c("#FF0000", "#7434a4"),
                     label = c(expression(paste(UP," ","in"," ",Foxp3^{Cre-YFP})),
                               expression(paste(UP," ","in"," ",Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
  # facet_wrap(~cellType, ncol = 3) + 
  ylab("Significance (-logPval)") +
  xlab("Log fold change") +
  # xlim(-2, 2) +
  theme_bw() +
  # guides(col=guide_legend(nrow=2,byrow=TRUE),
  #        text = NULL) +
  theme(text = element_text(family = "Arial", size = 18, face = "bold"),
        panel.border = element_rect(size = 2),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key=element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"))
dev.off()
