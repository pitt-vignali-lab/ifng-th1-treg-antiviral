# ----------------------------------------------------------------------------
# Fig 2a | UMAP of Treg clusters
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2A.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/2A.png
# Upstream : 02_treg_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(extrafont)
library(ggforce)

remotes::install_version("Rttf2pt1", version = "1.3.8")

font_import()

loadfonts(device = "win", quiet = TRUE)

fonts()

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

DimPlot(datLab, label = T) + NoLegend()

df <- data.frame(x = datLab[["umap"]]@cell.embeddings[,1],
                 y = datLab[["umap"]]@cell.embeddings[,2],
                 spl = datLab$samples,
                 lab = datLab@active.ident)

dfLab <- df %>% group_by(lab) %>% 
  summarise(x = mean(x), y = mean(y))

theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- cbind(cos(theta), sin(theta))
library(plyr)
aux <- function(x, one, two, prob=0.8) {
  if(nrow(x) <= 2) {
    return(NULL)
  }
  sigma <- var(cbind(x[,one], x[,two]))
  mu <- c(mean(x[,one]), mean(x[,two]))
  ed <- sqrt(qchisq(prob, df = 2))
  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
}

ell <- plyr::ddply(df, "lab", aux, one="x", two="y", prob = 0.2)



png("3.Output/treg/2A.png", width = 10.5, height = 10.5, units = "in", res = 300)

ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6, size = 3) +
  geom_mark_ellipse(data=ell, 
                    aes(x=X1, y=X2, label=lab, group=lab),
                    inherit.aes = F,
                    color = NA,
                    label.fill = "transparent",
                    label.fontsize = 25,
                    label.buffer = unit(25, 'mm'),
                    con.cap = 0) +
  theme_void()+ 
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  theme(legend.position="none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=2),
        axis.title = element_text(size = 25, face = "bold"),
        text = element_text(family = "Arial", size = 15))

dev.off()
