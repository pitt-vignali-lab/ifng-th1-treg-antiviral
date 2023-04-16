# ----------------------------------------------------------------------------
# Analysis | GSEA enrichment plots for the Treg comparison
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/4.2.gseaVis.R
# Input    : 3.Output/treg/
# Output   : 3.Output/treg/
# Upstream : 06_treg_gsea.R
# Note     : Customised clusterProfiler gseaplot / gseaplot2
# ----------------------------------------------------------------------------

library(tidyverse)
library(reshape)


for (i in dir("3.Output/treg/")[grep("gseaC5.*rds", dir("3.Output/treg/"))]) {
  assign(gsub(".rds", "", i),
         readRDS(paste0("3.Output/treg/", i)) %>% 
           filter(padj < 0.05))
}


gseaMNes <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), 
       list(gseaC5Overall[,c("pathway", "NES")],
            gseaC5Effector[,c("pathway", "NES")],
            gseaC5Central[,c("pathway", "NES")]))

names(gseaMNes) <- c("Pathway", "Overall", "Effector", "Central")

gseaMSz <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), 
                   list(gseaC5Overall[,c("pathway", "size")],
                        gseaC5Effector[,c("pathway", "size")],
                        gseaC5Central[,c("pathway", "size")]))

names(gseaMSz) <- c("Pathway", "Overall", "Effector", "Central")

commPath <- na.omit(gseaMNes)$Pathway
commPath <-  gsub("_", " ", commPath)
  
gseaMNes <- melt(gseaMNes, id = "Pathway")

gseaMSz <- melt(gseaMSz, id = "Pathway")

df <- cbind(gseaMNes, gseaMSz)
df <- df[c(1,2,3,6)]
names(df)[c(2,3,4)] <- c("Category", "NES", "Size")

df <- df[-grep("_UP|_DN", df$Pathway),]
df$Pathway <- gsub("_", " ", df$Pathway)
df$newPathway = str_wrap(df$Pathway, width = 45)
df$Category <- factor(df$Category, levels = c("Overall", "Central", "Effector"))

df <- df[grep("GOBP", df$Pathway),]
png("3.Output/treg/gseaC5CommBP.png", width = 8, height = 10, units = "in", res = 300)
ggplot(df[df$Pathway %in% commPath,], aes(x = Category, y = reorder(newPathway, NES), col = NES, size = Size)) +
  geom_point()+
  scale_color_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-4,4), space = "Lab") +
  ylab("Pathway") +
  facet_grid(~Category, scales="free_x")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

# enrichment plot
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
               "Metabolically_hyperactive", # 2
               "Proliferating", # 3
               "Gene_transcription_and_protein_stability" # 4
)

names(clusterId) <- levels(dat)

datLab <- RenameIdents(dat, clusterId)

Idents(datLab) <- "samples"

deOverall <- FindMarkers(datLab, ident.1 = "G", ident.2 = "C", 
                         test.use ="MAST",logfc.threshold = 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(p_val < 0.05)

library(msigdbr)
library(fgsea)
library(tibble)
library(viridis)


m_df_C2 <- msigdbr(species = "Mus musculus", category = "C2")
m_df_C5 <- msigdbr(species = "Mus musculus", category = "C5")
m_df_C7 <- msigdbr(species = "Mus musculus", category = "C7")

m_df <- rbind(m_df_C2, m_df_C5, m_df_C7)

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
selPath <- read.csv("TregSelectedGSEA.csv")
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

png("3.Output/treg/gseaRidge.png", width = 15, height = 10, units = "in", res = 300)
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
result <- as.data.frame(fgseaRes[fgseaRes$pathway %in% 
                                   selPath[selPath$Cat == "Ontology gene sets",]$Pathway,])

result$ID <- result$pathway
rownames(result) <- result$pathway
names(result)[c(1)] <- c("Description")

edo2m@result <- result

edo2m@geneSets <- fgsea_sets

edo2m@geneList <- ranks

# trace("gseaplot2", edit=TRUE)

png("3.Output/treg/gseaC5.png", width = 15, height = 8, units = "in", res = 300)
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
dir.create("3.Output/treg/seperateGSEA")

for (i in 1:nrow(edo2m)) {
  # png(paste0("3.Output/treg/seperateGSEA/gsea_", edo2m$Description[i], ".png"), width = 10, height = 6, units = "in", res = 300)
  p <- gseaplot2(edo2m, geneSetID = i) 
  ggsave(paste0("3.Output/treg/seperateGSEA/gsea_", edo2m$Description[i], ".png"), width = 8, height = 5, units = "in")
}

gsearank(edo2m, 10)
