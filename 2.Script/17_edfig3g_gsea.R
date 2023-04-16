# ----------------------------------------------------------------------------
# ED Fig 3g | GSEA of the Treg genotype comparison
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/S2D.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/S2D.png
# Upstream : 04_treg_de_mast.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(tibble)
library(extrafont)
library(ggnewscale)

dat <- readRDS("1.Data/treg/TregDat.RDS")
Idents(dat) <- "samples"
deOverall <- FindMarkers(dat, ident.1 = "G", ident.2 = "C", 
                         test.use ="MAST",logfc.threshold = 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(p_val < 0.05)

m_df<- msigdbr(species = "Mus musculus", category = "C2")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

gsea_genes<- deOverall %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(gsea_genes)

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=15,
                  maxSize=500)

esDf <- function(pathways, ids, stats, height = 0.5, gseaParam = 1) {
  df <- data.frame()
  for (i in 1:length(ids)) {
    pathway <- pathways[[ids[i]]]
    rnk <- rank(-ranks)
    ord <- order(rnk)
    statsAdj <- ranks[ord]
    
    statsAdj <- sign(statsAdj) * (abs(statsAdj))
    
    statsAdj <- statsAdj/max(abs(statsAdj))
    
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    
    toPlot$tops <- max(tops)
    toPlot$bottoms <- min(bottoms)
    toPlot$ID <- ids[i]
    
    df <- rbind(df, toPlot)
  }
  
  df$segmenty <- NA
  df$segmenty2 <- NA
  
  space <- sum(abs(range(df$tops, df$bottoms)))*height
  
  spaceAver <- seq(min(df$bottoms) *1.05, min(df$bottoms)*1.05 - space, -space/length(ids))
  
  for (i in length(ids):1) {
    df[df$ID == ids[i],]$segmenty = spaceAver[length(ids) - i + 1]
    df[df$ID == ids[i],]$segmenty2 = spaceAver[length(ids) - i + 2]
  }
  
  return(df)
}

selPath <- read.csv("TregSelectedGSEA.csv")$Pathway[11:13]

mydf <- esDf(fgsea_sets, selPath, ranks, height = 1/5)

mydf$ID <- sapply(mydf$ID, function(x) stringr::str_to_sentence(gsub("_", " ", x))) 
mydf$ID <- stringr::str_wrap(mydf$ID, 25)

names(mydf)[5] <- "Pathway"

range(ranks)

col_fun = circlize::colorRamp2(c(-0.9295100, 0, 0.6956676), c("#7434a4", "white", "#fb0707"))

png("3.Output/treg/S2D.png", width = 8.5*(0.75^0.5), height = 6*(0.75^0.5), units = "in", res = 300)
ggplot(mydf, aes(x = x, y = y, col = Pathway)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, 
             colour = "black",
             linetype = "dashed") +
  geom_segment(aes(x = x,y=segmenty,yend=segmenty2,xend=x)) +
  geom_rug(data=data.frame(x = (1:length(ranks)) - 1, y = as.numeric(ranks)),
           aes(x = x), sides="b",
           colour = ifelse(data.frame(x = (1:length(ranks)) - 1, y = as.numeric(ranks))$y > 0,
                           "#fb0707", "#7434a4"),
           length = unit(0.1, "npc")) +
  ylim(-1, 0.25) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 15, face = "bold")) +
  labs(x = "Rank", 
       y = "Enrichment score")
dev.off()

# ggplot(toPlot, aes(x = x, y = y)) + 
#   geom_line(color = "green")  + 
#   geom_rug(sides="b", col = "red")
#   geom_hline(yintercept = max(tops), colour = "red", 
#              linetype = "dashed") + 
#   geom_hline(yintercept = min(bottoms), 
#              colour = "red", linetype = "dashed") + 
#   geom_hline(yintercept = 0, 
#              colour = "black") + 
#   geom_line(color = "green") + 
#   theme_bw() + 
#   geom_segment(data = data.frame(x = pathway), 
#                mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), 
#                size = 0.2) + 
#   theme(panel.border = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   labs(x = "rank", 
#        y = "enrichment score")
