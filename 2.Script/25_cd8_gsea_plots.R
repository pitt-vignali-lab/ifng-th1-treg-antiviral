# ----------------------------------------------------------------------------
# Analysis | GSEA enrichment and ridge plots for CD8 T cells
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/4.2.gseaVis.R
# Input    : 3.Output/cd8/
# Output   : 3.Output/cd8/gseaC5.png
#            3.Output/cd8/gseaRidge.png
# Upstream : 24_cd8_gsea.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)

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
Idents(datLab) <- "samples"

deOverall <- FindMarkers(datLab, ident.1 = "G", ident.2 = "C", 
                         test.use ="MAST",logfc.threshold = 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(p_val < 0.05)

library(msigdbr)
library(fgsea)
library(tibble)
library(viridis)
library(DOSE)
library(enrichplot)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo2 <- gseDO(geneList)

# m_df_C2 <- msigdbr(species = "Mus musculus", category = "C2")
m_df_C5 <- msigdbr(species = "Mus musculus", category = "C5")
# m_df_C7 <- msigdbr(species = "Mus musculus", category = "C7")

m_df <- rbind(#m_df_C2, 
              m_df_C5)#, 
              # m_df_C7)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

gsea_genes<- deOverall %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(gsea_genes)

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=15,
                  maxSize=500)

# ridges ----
selPath <- read.csv("selectedGSEA.csv")
gs2id <- fgsea_sets[selPath$Pathway]

gs2val <- lapply(gs2id, function(id) {
  res <- ranks[id]
  res <- res[!is.na(res)]
})

nn <- names(gs2val)
i <- match(nn, fgseaRes$pathway)

j <- order(fgseaRes[["padj"]][i], decreasing = F)
len <- sapply(gs2val, length)
gs2val.df <- data.frame(category = rep(nn, times=len),
                        color = rep(as.data.frame(fgseaRes)[i,"padj"], times=len),
                        value = unlist(gs2val))

colnames(gs2val.df)[2] <- "padj"

gs2val.df <- merge(selPath, gs2val.df, by.x = "Pathway", by.y = "category")

png("3.Output/cd8/gseaRidge.png", width = 15, height = 5, units = "in", res = 300)
ggplot(gs2val.df, aes(x=value, y=reorder(Pathway, padj), fill=padj)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  scale_fill_viridis(name = "Adjusted P value",
                     guide=guide_colorbar(reverse=T)) +
  # scale_fill_continuous(low="yellow", high="orange", name = "padj", 
  #                       guide=guide_colorbar(reverse=T)) +
  facet_grid(Cat~., switch = "y", scales = "free_y", space = "free_y") +
  xlab("Gene enrichment") + 
  ylab(NULL)  +
  theme_bw() +
  theme(text = element_text(size=17, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, 'cm'),
        strip.background = element_rect(colour = "black", fill = "white"))
dev.off()

### combined gsea plot---
edo2m <- edo2
result <- as.data.frame(fgseaRes[fgseaRes$pathway %in% 
                                   selPath[selPath$Cat == "Ontology gene sets",]$Pathway,])

result$ID <- result$pathway
rownames(result) <- result$pathway
names(result)[c(1)] <- c("Description")

edo2m@result <- result

edo2m@geneSets <- fgsea_sets

edo2m@geneList <- ranks

# trace("gseaplot2", edit=TRUE)

png("3.Output/cd8/gseaC5.png", width = 15, height = 8, units = "in", res = 300)
gseaplot2(edo2m, geneSetID = 1:nrow(result)) 
dev.off()

### seperate gsea plot---

result <- as.data.frame(fgseaRes[fgseaRes$pathway %in% selPath$Pathway,])

result$ID <- result$pathway
rownames(result) <- result$pathway
names(result)[c(1)] <- c("Description")

edo2m@result <- result

edo2m@geneSets <- fgsea_sets

edo2m@geneList <- ranks

# trace("gseaplot2", edit=TRUE)
dir.create("3.Output/cd8/seperateGSEA")

for (i in 1:nrow(edo2m)) {
  # png(paste0("3.Output/cd8/seperateGSEA/gsea_", edo2m$Description[i], ".png"), width = 10, height = 6, units = "in", res = 300)
  p <- gseaplot2(edo2m, geneSetID = i) 
  ggsave(paste0("3.Output/cd8/seperateGSEA/gsea_", edo2m$Description[i], ".png"), width = 8, height = 5, units = "in")
}

gsearank(edo2m, 10)
