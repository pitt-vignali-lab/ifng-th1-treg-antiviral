# ----------------------------------------------------------------------------
# Revision | GSEA and SPIA on the pseudobulk comparison
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/5.3.gsea&spia.R
# Input    : 3.Output/revision/4.PseudoBulk/
# Output   : 3.Output/revision/5.PseudoBulk/SPIA_KEGG.csv
# Upstream : 35_rev_pseudobulk.R
# Note     : SPIA needs a local KEGG data directory at 1.Data/revision/RefData/
# ----------------------------------------------------------------------------

library(msigdbr)
library(fgsea)
library(tibble)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(SPIA)

rank <- read.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv",
               row.names = 1) %>% 
  arrange(desc(logFC)) %>%
  filter(logFC != 0) %>% 
  dplyr::select(feature, logFC) %>% 
  deframe()

# GSEA ----
fgsea_sets<- msigdbr(species = "Mus musculus", 
                     category = "C5", # C2, C7
                     subcategory = "BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name) 

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = rank,
                  minSize=5,
                  nproc = 125)

fgseaRes[grep("_T_CELL_MIG", fgseaRes$pathway),]

fgseaRes %>% filter(padj < 0.05) %>% nrow()

data.table::fwrite(fgseaRes, "3.Output/revision/5.PseudoBulk/GSEA_GOBP_G.csv")

# SPIA ----

de <- read.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv",
               row.names = 1)

allGene <- bitr(de$feature, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db)

genes <- merge(de, allGene, by.x = "feature", by.y = "SYMBOL")

ranks <- genes %>% 
  filter(padj < 0.05) %>% 
  dplyr::select( "ENTREZID", "logFC") %>% 
  arrange(desc(logFC)) %>% 
  deframe()

spiaRes <- spia(de = ranks, 
                all = allGene$ENTREZID, 
                organism = "mmu", 
                data.dir = "1.Data/revision/RefData/", 
                nB = 2000)
 spiaRes %>% filter(pGFdr < 0.05) %>% nrow()

write.csv(spiaRes, "3.Output/revision/5.PseudoBulk/SPIA_KEGG.csv") 

