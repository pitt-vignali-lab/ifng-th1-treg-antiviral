# ----------------------------------------------------------------------------
# Revision | GSEA enrichment plots for the pseudobulk comparison
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/5.4.gseaPlot.R
# Input    : 3.Output/revision/5.PseudoBulk/
# Output   : 3.Output/revision/5.PseudoBulk/
# Upstream : 37_rev_gsea.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(tibble)
library(extrafont)
library(ggnewscale)


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

resDf <- read.csv("3.Output/revision/5.PseudoBulk/GSEA_IMMUNE_SIGNAT_G.csv")

# individual plots ----
ranks <- read.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv",
                 row.names = 1) %>% 
  arrange(desc(logFC)) %>%
  filter(logFC != 0) %>% 
  dplyr::select(feature, logFC) %>% 
  deframe()

m_df<- msigdbr(species = "Mus musculus", category = "C7")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks,
                  minSize=15,
                  maxSize=500)

for( selPath in c("GSE39864_WT_VS_GATA3_KO_TREG_UP", #3C ylim(-0.8, 0)
                  "GSE17974_CTRL_VS_ACT_IL4_AND_ANTI_IL12_72H_CD4_TCELL_UP", #3B ylim(-0.2, 0.5)
                  "GSE40685_TREG_VS_FOXP3_KO_TREG_PRECURSOR_UP",
                  "GSE34006_A2AR_KO_VS_A2AR_AGONIST_TREATED_TREG_UP",
                  "GSE17974_IL4_AND_ANTI_IL12_VS_UNTREATED_24H_ACT_CD4_TCELL_DN",
                  "GSE32423_CTRL_VS_IL4_MEMORY_CD8_TCELL_DN",
                  "GSE40443_INDUCED_VS_TOTAL_TREG_DN",
                  "GSE14415_INDUCED_TREG_VS_FOXP3_KO_INDUCED_TREG_UP",
                  "GSE40685_TREG_VS_FOXP3_KO_TREG_PRECURSOR_UP",
                  "GSE7852_TREG_VS_TCONV_UP")){
  
  mydf <- esDf(fgsea_sets, selPath, ranks, height = 1/7)
  
  mydf$ID <- sapply(mydf$ID, function(x) gsub("_", " ", x))
  
  names(mydf)[5] <- "Pathway"
  
  range(ranks)
  
  ggplot(mydf, aes(x = x, y = y)) +
    geom_line(linewidth = 1.5, color = "green") +
    geom_hline(yintercept = 0, 
               colour = "black",
               linetype = "dashed") +
    geom_segment(aes(x = x, y=segmenty, 
                     yend=segmenty2, xend=x,
                     col = ifelse(x < 8845, "UP", "DN")),
                 linewidth = 0.7, alpha = 0.3) +
    scale_color_manual(breaks = c("UP", "DN"),
                       values = c("#7434a4", "#fb0707"),
                       labels = expression(bold(Ifngr1^{L/L}~Foxp3^{Cre/YFP}), 
                                           bold(Foxp3^{Cre/YFP}))) +
    guides(color = guide_legend(
      override.aes = list(alpha = 1, 
                          linewidth = 7))) +
    labs(col = "Genotype") +
    geom_vline(xintercept = 8845,
               colour = "gray80",
               linetype = "dashed") +
    theme_bw()  +
    labs(x = "Rank", 
         y = "Enrichment score") +
    ggtitle(mydf$Pathway, 
            subtitle = paste("NES:", 
                             format(fgseaRes[fgseaRes$pathway == selPath,]$NES, 
                                    digits = 3),
                             "padj:", 
                             format(fgseaRes[fgseaRes$pathway == selPath,]$padj, 
                                    digits = 3, scientific = T))) +
    theme(legend.position = "bottom",
          text = element_text(family = "sans", size = 15, face = "bold"),
          panel.border = element_rect(colour = "black", linewidth = 1.5),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10, face = "bold"))
  
  ggsave(paste0("3.Output/revision/5.PseudoBulk/", 
                strsplit(selPath, split = "_")[[1]][1],
                ".png"), 
         width = 6, height = 4,
         units = "in", bg = "white")
  
} #4B ylim(-0.9, 0)


# combined plots ----
selPath <- c("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
             "BOSCO_TH1_CYTOTOXIC_MODULE",
             "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM")

m_df<- msigdbr(species = "Mus musculus", category = "C2") %>% 
  filter(gs_name %in% selPath)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks,
                  minSize=15,
                  maxSize=1000)

mydf <- esDf(fgsea_sets, selPath, ranks, height = 1/4)

mydf$ID <- sapply(mydf$ID, function(x) stringr::str_to_sentence(gsub("_", " ", x))) 

annoDf <- data.frame(ID = stringr::str_to_sentence(gsub("_", " ", fgseaRes$pathway)),
                     ID2 = paste0(stringr::str_to_sentence(gsub("_", " ", fgseaRes$pathway)),
                                  " NES:",
                                  format(fgseaRes$NES, digits = 3),
                                  " padj:",
                                  format(fgseaRes$padj, digits = 3, scientific = T)))

mydf <- merge(mydf, annoDf, by = "ID")
# mydf$ID <- stringr::str_wrap(mydf$ID, 25)

names(mydf)[5] <- "Pathway"

range(ranks)

col_fun = circlize::colorRamp2(c(-0.9295100, 0, 0.6956676), c("#7434a4", "white", "#fb0707"))

ggplot(mydf, aes(x = x, y = y, col = ID2)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, 
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = 8845,
             colour = "gray80",
             linetype = "dashed") +
  geom_segment(aes(x = x,y=segmenty,yend=segmenty2,xend=x)) +
  geom_rug(data=data.frame(x = (1:length(ranks)) - 1, y = as.numeric(ranks)),
           aes(x = x), sides="b",
           colour = ifelse(data.frame(x = (1:length(ranks)) - 1, y = as.numeric(ranks))$y > 0,
                           "#7434a4", "#fb0707"),
           length = unit(0.03, "npc")) +
  ylim(-1.15, 0.25) +
  theme_bw() +
  guides(color = guide_legend(
    nrow=3, byrow=TRUE,
    override.aes = list(alpha = 1, 
                        linewidth = 7))) +
  labs(col = "",) +
  theme(legend.position = "bottom",
        text = element_text(family = "sans", size = 12.5, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 1.5),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "bold")) +
  labs(x = "Rank", 
       y = "Enrichment score")

ggsave("3.Output/revision/5.PseudoBulk/GSEAcombined.png", 
       width = 6, height = 4.5,
       units = "in", bg = "white")

## flip over version ----
ranks <- read.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv",
                  row.names = 1) %>% 
  arrange(desc(logFC)) %>%
  filter(logFC != 0) %>% 
  dplyr::select(feature, logFC) %>% 
  deframe()

ranks <- ranks * -1

m_df<- msigdbr(species = "Mus musculus", category = "C7")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks,
                  minSize=15,
                  maxSize=500)

for( selPath in c("GSE39864_WT_VS_GATA3_KO_TREG_UP", #3C ylim(-0.8, 0)
                  "GSE17974_CTRL_VS_ACT_IL4_AND_ANTI_IL12_72H_CD4_TCELL_UP", #3B ylim(-0.2, 0.5)
                  "GSE40685_TREG_VS_FOXP3_KO_TREG_PRECURSOR_UP",
                  "GSE34006_A2AR_KO_VS_A2AR_AGONIST_TREATED_TREG_UP",
                  "GSE17974_IL4_AND_ANTI_IL12_VS_UNTREATED_24H_ACT_CD4_TCELL_DN",
                  "GSE32423_CTRL_VS_IL4_MEMORY_CD8_TCELL_DN",
                  "GSE40443_INDUCED_VS_TOTAL_TREG_DN",
                  "GSE14415_INDUCED_TREG_VS_FOXP3_KO_INDUCED_TREG_UP",
                  "GSE40685_TREG_VS_FOXP3_KO_TREG_PRECURSOR_UP",
                  "GSE7852_TREG_VS_TCONV_UP")){
  
  mydf <- esDf(fgsea_sets, selPath, ranks, height = 1/7)
  
  mydf$ID <- sapply(mydf$ID, function(x) gsub("_", " ", x))
  
  names(mydf)[5] <- "Pathway"
  
  range(ranks)
  
  ggplot(mydf, aes(x = x, y = y)) +
    geom_line(linewidth = 1.5, color = "green") +
    geom_hline(yintercept = 0, 
               colour = "black",
               linetype = "dashed") +
    geom_segment(aes(x = x, y=segmenty, 
                     yend=segmenty2, xend=x,
                     col = ifelse(x < 9226, "UP", "DN")),
                 linewidth = 0.7, alpha = 0.3) +
    scale_color_manual(breaks = c("UP", "DN"),
                       values = c("#fb0707", "#7434a4"),
                       labels = expression(bold(Foxp3^{Cre/YFP}),
                                           bold(Ifngr1^{L/L}~Foxp3^{Cre/YFP}))) +
    guides(color = guide_legend(
      override.aes = list(alpha = 1, 
                          linewidth = 7))) +
    labs(col = "Genotype") +
    geom_vline(xintercept = 9226 - 1,
               colour = "gray80",
               linetype = "dashed") +
    theme_bw()  +
    labs(x = "Rank", 
         y = "Enrichment score") +
    ggtitle(mydf$Pathway, 
            subtitle = paste("NES:", 
                             format(resDf[resDf$pathway == selPath,]$NES * -1, 
                                    digits = 3),
                             "padj:", 
                             format(resDf[resDf$pathway == selPath,]$padj, 
                                    digits = 3, scientific = T))) +
    theme(legend.position = "bottom",
          text = element_text(family = "sans", size = 15, face = "bold"),
          panel.border = element_rect(colour = "black", linewidth = 1.5),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10, face = "bold"))
  
  ggsave(paste0("3.Output/revision/5.PseudoBulk/flipOver_", 
                selPath,
                ".png"), 
         width = 6, height = 4,
         units = "in", bg = "white")
  
}

########## combined #################
selPath <- c("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
             "BOSCO_TH1_CYTOTOXIC_MODULE",
             "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM")

m_df<- msigdbr(species = "Mus musculus", category = "C2") %>% 
  filter(gs_name %in% selPath)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks,
                  minSize=15,
                  maxSize=1000)

mydf <- esDf(fgsea_sets, selPath, ranks, height = 1/4)

mydf$ID <- sapply(mydf$ID, function(x) stringr::str_to_sentence(gsub("_", " ", x))) 

annoDf <- data.frame(ID = stringr::str_to_sentence(gsub("_", " ", fgseaRes$pathway)),
                     ID2 = paste0(stringr::str_to_sentence(gsub("_", " ", fgseaRes$pathway)),
                                  " NES:",
                                  format(fgseaRes$NES, digits = 3),
                                  " padj:",
                                  format(fgseaRes$padj, digits = 3, scientific = T)))

mydf <- merge(mydf, annoDf, by = "ID")
# mydf$ID <- stringr::str_wrap(mydf$ID, 25)

names(mydf)[5] <- "Pathway"

range(ranks)

col_fun = circlize::colorRamp2(c(-0.9295100, 0, 0.6956676), c("#7434a4", "white", "#fb0707"))

ggplot(mydf, aes(x = x, y = y, col = ID2)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, 
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = table(ranks > 0)[2] - 1,
             colour = "gray80",
             linetype = "dashed") +
  geom_segment(aes(x = x,y=segmenty,yend=segmenty2,xend=x)) +
  geom_rug(data=data.frame(x = (1:length(ranks)) - 1, 
                           y = sort(as.numeric(ranks),
                                    decreasing = T)),
           aes(x = x), sides="b",
           colour = ifelse(data.frame(x = (1:length(ranks)) - 1, 
                                      y = sort(as.numeric(ranks),
                                               decreasing = T))$y > 0,
                           "#fb0707", "#7434a4"),
           length = unit(0.03, "npc")) +
  ylim(-0.5, 0.9) +
  theme_bw() +
  guides(color = guide_legend(
    nrow=3, byrow=TRUE,
    override.aes = list(alpha = 1, 
                        linewidth = 7))) +
  labs(col = "",) +
  theme(legend.position = "bottom",
        text = element_text(family = "sans", size = 12.5, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 1.5),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "bold")) +
  labs(x = "Rank", 
       y = "Enrichment score")

ggsave("3.Output/revision/5.PseudoBulk/flipOver_GSEAcombined.png", 
       width = 6, height = 4.5,
       units = "in", bg = "white")
