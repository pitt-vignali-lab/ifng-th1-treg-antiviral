# ----------------------------------------------------------------------------
# Fig 2b | Distribution of control and Ifngr1-deficient Tregs across clusters
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2B.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/2B_Galaxy*.png
# Upstream : 02_treg_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(extrafont)
library(ggdark)
library(viridis)
library(ggforce)

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

get_density <- function(x, y, ...) {
  density_out <- MASS::kde2d(x, y, ...)
  int_x <- findInterval(x, density_out$x)
  int_y <- findInterval(y, density_out$y)
  comb_int <- cbind(int_x, int_y)
  return(density_out$z[comb_int])
}

df <- df %>% group_by(lab) %>% 
  mutate(density = get_density(x, y, h = c(1, 1), n = 100))


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

ell <- plyr::ddply(df, "lab", aux, one="x", two="y", prob = 0.3)


png("3.Output/treg/2B_GalaxyDark1.png", width = 15, height = 15, units = "in", res = 300)
ggplot(df) + aes(x=x, y=y) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
  geom_point(size=0.1) +
  stat_density2d(geom="tile", aes(fill=..density..^0.25,     
                                  alpha=ifelse(..density..^0.25<0.1,1,0)), contour=FALSE) + 
  scale_fill_viridis(option="magma") +
  geom_mark_ellipse(data=ell, aes(x=X1, y=X2, label=lab, col=lab),
                    inherit.aes = F, size = 2, 
                    label.fontsize = 19,
                    label.buffer = unit(2, 'mm')) +
  xlim(-10, 10) +
  ylim(-8, 8) +
  dark_theme_void() +
  facet_wrap(~spl, ncol = 2) +
  theme(legend.position = "none")
dev.off()

png("3.Output/treg/2B_GalaxyLight.png", width = 10, height = 10, units = "in", res = 300)
ggplot(df) + aes(x=x, y=y) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
  geom_point(size=0.1) +
  stat_density2d(geom="tile", aes(fill=..density..^0.25,     
                                  alpha=ifelse(..density..^0.25<0.1,1,0)), 
                 contour=FALSE) + 
  # scale_fill_viridis(option="magma") +
  xlim(-10, 10) +
  ylim(-8, 8) +
  theme_void() +
  facet_wrap(~spl, ncol = 2) +
  theme(legend.position = "none")
dev.off()