# ----------------------------------------------------------------------------
# Data mining | Cytokine-receptor pair interactions across LCMV-infected cell types
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/6.3.3.gse186898celltalker.R
# Input    : 1.Data/revision/mouseLRModified.csv
#            1.Data/revision/mouseChemokine*.csv
# Output   : 3.Output/revision/
# Upstream : 40_gse186898_cluster_annotation.R
# Note     : celltalker
#            pairs merged and deduplicated from several sources, P < 0.05 retained
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor")){
  dir.create("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor")
}

library(celltalker)
library(Seurat)
library(dplyr)
library(tidyverse)
library(nichenetr)

# reference data read in ----
refDf <- cbind(read.csv("1.Data/revision/mouseLRModified.csv", row.names = 1)%>% 
                 rename(from = 1, toHsa = 2),
               read.csv("1.Data/revision/mouseLRModified.csv", row.names = 1) %>% 
                 rename(from = 1, to = 2) %>% 
                 mutate(to = convert_human_to_mouse_symbols(to)) %>% 
                 select(to))

refDf %>% filter(is.na(to))

refDf[refDf$toHsa == "D6",]$to <- "Ackr2"
refDf[refDf$toHsa == "DARC",]$to <- "Ackr1"
refDf[refDf$toHsa == "CCR11",]$to <- "Ackr4" # Cclr2
refDf[refDf$toHsa == "CXCR7",]$to <- "Ackr3"
refDf[refDf$toHsa == "CXCR3B",]$to <- "Cxcr3"
refDf[refDf$toHsa == "FPRL-1",]$to <- "Fpr2"

refDf[refDf$from == "Cxcl4",]$from <- "Pf4"
refDf[refDf$from == "Ccl26l",]$from <- "Ccl26"
refDf[refDf$from == "Gm1960",]$from <- "Cxcl3"

refDf <- refDf %>% 
  unique() %>% 
  select(from, to) %>% 
  rename(ligand = 1, receptor = 2) %>% 
  mutate(pair = paste(ligand, receptor, sep = "_"))

# target data ----
dat <- readRDS("3.Output/revision/6.Chemokins/GSE186898/2.datAnnotated.RDS")

dat@meta.data <- dat@meta.data %>% 
  separate(title, sep = "_", into = c("genotype", "time"), remove = F)

refDf$ligand[!refDf$ligand %in% rownames(dat)]

# analysis running ----

int <- celltalk(input_object = dat,
                metadata_grouping = "anno",
                ligand_receptor_pairs = refDf,
                number_cells_required = 50,
                min_expression = 100,
                max_expression = 9000,
                scramble_times = 10)

intSub <- int %>% filter(cell_type1 == "Treg" | cell_type2 == "Treg") %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>% 
  filter(p_val<0.1) %>% 
  separate(interaction, sep = "_", into = c("ligand", "receptpor"), remove = F)

intSub %>% filter(fdr < 0.05) %>% nrow()

write.csv(intSub[c(9, 10, 3, 4, 8, 12, 11)], 
          "3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/Overall_lr.csv")

write.csv(refDf, "3.Output/revision/6.Chemokins/chemokinReceptorPairs.csv")

for (tm in unique(dat$time)) {
  int <- celltalk(input_object = subset(dat, time == tm),
                  metadata_grouping = "anno",
                  ligand_receptor_pairs = refDf,
                  number_cells_required = 50,
                  min_expression = 100,
                  max_expression = 9000,
                  scramble_times = 10)
  
  intSub <- int %>% filter(cell_type1 == "Treg" | cell_type2 == "Treg") %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>% 
    filter(fdr<0.05) %>% 
    separate(interaction, sep = "_", into = c("ligand", "receptpor"), remove = F)
  
  write.csv(intSub[c(9, 10, 3, 4, 8, 12, 11)], 
            paste0("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/", tm, "_lr.csv"))
}

for (tm in unique(dat$time)) {
  for (gt in unique(dat$genotype)) {
    int <- celltalk(input_object = subset(dat, time == tm & genotype == gt),
                    metadata_grouping = "anno",
                    ligand_receptor_pairs = lrSub,
                    number_cells_required = 50,
                    min_expression = 100,
                    max_expression = 20000,
                    scramble_times = 10)
    
    intSub <- int %>% filter(cell_type1 == "Treg" | cell_type2 == "Treg") %>%
      mutate(fdr=p.adjust(p_val,method="fdr")) %>% 
      filter(fdr<0.05) %>% 
      separate(interaction, sep = "_", into = c("ligand", "receptpor"), remove = F)
    
    write.csv(intSub[c(9, 10, 3, 4, 8, 12, 11)], 
              paste0("3.Output/revision/6.Chemokins/GSE186898/3.ligandReceptor/", tm, "_", gt, "_lr.csv")) 
  }
}

# visulize ----

colors_use <- RColorBrewer::brewer.pal(n=length(unique(c(intSub$cell_type1, intSub$cell_type2))),"Set2")

circos_plot(ligand_receptor_frame=intSub,
            cell_group_colors=colors_use,
            ligand_color="blue",
            receptor_color="red",
            cex_outer=0.5,
            cex_inner=0.4)
