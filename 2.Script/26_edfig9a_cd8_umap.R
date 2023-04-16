# ----------------------------------------------------------------------------
# ED Fig 9a | UMAP of tetramer-positive CD8 T cell clusters
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/7A.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/7A.png
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(extrafont)
library(ggforce)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")
DefaultAssay(dat) <- "RNA"
Idents(dat) <- "RNA_snn_res.0.3"

datSub <- subset(dat, idents = 8, invert = TRUE)

clusterId <- c("Activated", # 0
               "Cytotoxic", # 1
               "Progenitor\n exhausted", # 2
               "Proliferating", # 3
               "Exhausted", # 4
               "Memory", # 5
               "IFN signaling", # 6
               "Metabolically\n hyperactive" # 7
)

names(clusterId) <- levels(datSub)

datLab <- RenameIdents(datSub, clusterId)

DimPlot(datLab, label = T) + NoLegend()


df <- data.frame(x = datLab[["umap"]]@cell.embeddings[,1],
                 y = datLab[["umap"]]@cell.embeddings[,2],
                 lab = datLab@active.ident,
                 rep = datLab$reps)

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

png("3.Output/cd8/7A.png", width = 10, height = 10, units = "in", res = 300)

ggplot(df, aes(x , y, col = lab, fill = lab)) +
  geom_point(shape=21, alpha = 0.6, size = 3) +
  geom_mark_ellipse(data=ell, 
                    aes(x=X1, y=X2, label=lab, group=lab),
                    inherit.aes = F,
                    color = NA,
                    label.fill = "transparent",
                    label.fontsize = 25,
                    # con.border = 'none', 
                    # con.type = "none",
                    label.buffer = unit(3, 'mm'),
                    con.cap = 0) +
  theme_void()+ 
  # xlim(-8, 8)+
  # ylim(-6.5, 6.5) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  theme(legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 25, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(family = "Arial", size = 15))

dev.off()
