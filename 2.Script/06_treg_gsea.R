# ----------------------------------------------------------------------------
# Analysis | GSEA and SPIA pathway analysis of the Treg DE results
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/4.1.gseaSpia.R
# Input    : 3.Output/treg/deMAST_*.csv
# Output   : 3.Output/treg/SPIA_*.csv
# Upstream : 04_treg_de_mast.R
# Note     : MSigDB C2 and C7 gene sets, set size 15-500
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
               "Metabolically_hyperactive", # 2
               "Proliferating", # 3
               "Gene_transcription_and_protein_stability" # 4
)

names(clusterId) <- levels(dat)

datLab <- RenameIdents(dat, clusterId)

DimPlot(datLab, label = T) + NoLegend()

# 3. GSEA ----
library(msigdbr)
library(fgsea)
library(tibble)


# de ----
datLab <- RenameIdents(dat, clusterId)

for (i in clusterId) {
  datSub <- subset(datLab, idents = i)
  Idents(datSub) <- "samples"
  assign(paste0("de", i), 
         FindMarkers(datSub, ident.1 = "G", ident.2 = "C",
                     test.use = "MAST", logfc.threshold =  0.05) %>%
           rownames_to_column(var = "gene") %>% 
           filter(p_val < 0.05))
}

Idents(datLab) <- "samples"

deOverall <- FindMarkers(datLab, ident.1 = "G", ident.2 = "C", 
                         test.use ="MAST",logfc.threshold = 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(p_val < 0.05)

# GSEA 
cat ="C7"

m_df<- msigdbr(species = "Mus musculus", category = cat)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

for (j in ls()[grep("de", ls())]) {
  gsea_genes<- eval(parse(text = j)) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(gsea_genes)
  
  fgseaRes <- fgsea(pathways = fgsea_sets,
                    stats = ranks ,
                    minSize=15,
                    maxSize=500)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    filter(padj < 0.05)
  
  saveRDS(fgseaRes, paste0("3.Output/treg/gsea", cat, gsub("de","",j), ".rds"))
  write.csv(apply(fgseaResTidy,2,as.character), 
            paste0("3.Output/treg/gsea", cat, gsub("de","",j), ".csv"))
}

# SPIA ----
library(clusterProfiler)
library(org.Mm.eg.db)
library(SPIA)

allGene <- bitr(rownames(datLab), fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db)$ENTREZID

for (j in ls()[grep("de", ls())][-3]) {
  genes <- eval(parse(text = j))
  
  gRef <- bitr(genes$gene, fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  
  genes <- merge(genes, gRef, by.x = "gene", by.y= "SYMBOL")
  ranks <- genes[genes$p_val_adj<0.05,]$avg_log2FC
  names(ranks) <- genes[genes$p_val_adj<0.05,]$ENTREZID
  ranks <- sort(ranks, decreasing = T)

  # allGene <- unique(genes$ENTREZID)
  spiaRes <- spia(de = ranks, all = allGene, organism = "mmu", data.dir = "./", nB = 2000)
  
  write.csv(spiaRes, 
            paste0("3.Output/treg/spia", gsub("de","",j), ".csv"))}

# fgseaResTidy %>%
#   dplyr::select(-leadingEdge, -ES) %>%
#   arrange(padj) %>%
#   head()
# 
# ggplot(fgseaResTidy %>% filter(pval < 0.5) %>% head(n= 15), aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES < 0)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA") +
#   theme_minimal() 
# 
# plotEnrichment(fgsea_sets[["GRAESSMANN_APOPTOSIS_BY_DOXORUBICIN_UP"]],
#                ranks) + labs(title="GRAESSMANN_APOPTOSIS_BY_DOXORUBICIN_UP")
