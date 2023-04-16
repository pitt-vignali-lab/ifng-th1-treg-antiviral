# ----------------------------------------------------------------------------
# Analysis | GSEA and SPIA pathway analysis of the CD8 DE results
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/4.1.gseaSpia.R
# Input    : 3.Output/cd8/deMAST_Overall.csv
# Output   : 3.Output/cd8/
# Upstream : 23_cd8_de_analysis.R
# ----------------------------------------------------------------------------

# 1. Data initialization and labelling ----
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
cat ="C2"

m_df<- msigdbr(species = "Mus musculus", category = cat)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

for (j in ls()[grep("de", ls())]) {
  gsea_genes<- eval(parse(text = paste0("`", j, "`"))) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  
  assign(paste0(gsub("de", "",j), "_ranks"), deframe(gsea_genes))
  
  assign(paste0(gsub("de", "",j), "_fgseaRes"), 
         fgsea(pathways = fgsea_sets,
               stats = eval(parse(text = paste0("`", gsub("de", "",j), "_ranks`"))),
               minSize=15,
               maxSize=500))
  
  fgseaResTidy <- eval(parse(text = paste0("`", gsub("de", "",j), "_fgseaRes`"))) %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    filter(padj < 0.05)
  
  save(list = c(paste0(gsub("de", "",j), "_ranks"),
                paste0(gsub("de", "",j), "_fgseaRes")), 
          file = paste0("3.Output/cd8/gsea", cat, gsub("de","",j), ".RData"))
  write.csv(apply(fgseaResTidy,2,as.character), 
            paste0("3.Output/cd8/gsea", cat, gsub("de","",j), ".csv"))
}

plotEnrichment(
  fgsea_sets[["OSMAN_BLOOD_CHAD63_KH_AGE_18_50YO_HIGH_DOSE_SUBJECTS_24HR_UP"]], deActivated_ranks) 

# SPIA ----
library(clusterProfiler)
library(org.Mm.eg.db)
library(SPIA)

allGene <- bitr(rownames(datLab), fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db)$ENTREZID

for (j in ls()[grep("de", ls())][9]) {
  genes <- eval(parse(text = paste0("`", j, "`")))
  
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
            paste0("3.Output/cd8/spia", gsub("de","",j), ".csv"))}

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
