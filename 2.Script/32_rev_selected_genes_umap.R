# ----------------------------------------------------------------------------
# Revision | Selected gene expression on the integrated UMAP
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.2.comparesionSelectedGenesInUmap.R
# Input    : 1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS
# Output   : 3.Output/revision/2.Integration/
# Upstream : 31_rev_batch_correction_simspec.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggpubr)
library(presto)

dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")
DimPlot(dat, split.by = "samples") + coord_fixed()
# marker genes ineach cluster ----
deDf <- wilcoxauc(dat, group_by = "newAnno")

deDf %>% 
  filter(logFC > 0) %>% 
  group_by(group) %>%
  top_n(n = 20, wt = -padj) %>% 
  write.csv("3.Output/revision/2.Integration/top20MarkGenesByPadjEachCluster.csv",
            row.names = F)

deDf %>% 
  filter(padj < 0.05) %>% 
  group_by(group) %>%
  top_n(n = 20, wt = logFC) %>% 
  write.csv("3.Output/revision/2.Integration/top20MarkGenesByLogfcEachCluster.csv",
            row.names = F)

# combination umaps ----
cDat <- subset(dat, samples == "C")
gDat <- subset(dat, samples == "G")

p <- list()

for (i in strsplit("Stat1, Cxcr3, Ccl5, Il16, Nkg7, Klrg1, Ifng, Nfkbia, Il1rl1, Il7r, Gata3, Tbx21, Gpx4", ", ")[[1]]) {
  
  pc <- FeaturePlot(cDat, features = i, combine = FALSE, coord.fixed = T) 
  pc <- lapply(pc, function(x) x + 
                 scale_color_gradientn(colours = c("gray90", "#F95700FF"),  limits = c(0, 5))+
                 labs(subtitle = "C1-5") +
                 coord_fixed() +
                 theme(plot.subtitle = element_text(hjust = 0.5)))
  
  pg <- FeaturePlot(gDat, features = i, combine = FALSE, coord.fixed = T) 
  pg <- lapply(pg, function(x) x + 
                 scale_color_gradientn(colours = c("gray90", "#4203c9"),  limits = c(0, 5))+
                 labs(subtitle = "G1-5") +
                 coord_fixed() +
                 theme(plot.subtitle = element_text(hjust = 0.5)))
  
  p[[i]] <- CombinePlots(pc) + CombinePlots(pg)
}

png("3.Output/revision/combinedGenes.png", width = 20, height = 15, units = "in", res = 300)
CombinePlots(p, ncol = 3)
dev.off()
