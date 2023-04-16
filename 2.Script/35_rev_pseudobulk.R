# ----------------------------------------------------------------------------
# Revision | Pseudobulk aggregation by sample and cluster
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/5.1.pseudoBulk.R
# Input    : 1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS
# Output   : 3.Output/revision/4.PseudoBulk/
# Upstream : 31_rev_batch_correction_simspec.R
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/5.PseudoBulk/")){
  dir.create("3.Output/revision/5.PseudoBulk/")
}

library(Seurat)
library(DESeq2)
library(sva)
library(tidyverse)
library(presto)

dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")

wilcoxauc(dat, group_by = "samples", 
          groups_use = c("C", "G")) %>%
  filter(group == "G") %>%
  write.csv("3.Output/revision/5.PseudoBulk/deDfGvsC.csv")

wilcoxauc(dat, group_by = "samples", 
          groups_use = c("UC", "UG")) %>%
  filter(group == "UG") %>%
  write.csv("3.Output/revision/5.PseudoBulk/deDfUGvsUC.csv")


