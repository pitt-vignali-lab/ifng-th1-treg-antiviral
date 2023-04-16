# ----------------------------------------------------------------------------
# ED Fig 9a | Assign CD8 T cell cluster identities from marker genes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/2.1.cellTypeLabeling.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/2.top20SigMarkerRes0.3.csv
# Upstream : 19_cd8_data_init.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")
Idents(dat) <- "RNA_snn_res.0.3"
DimPlot(dat, label = T) +
  NoLegend()

markers <- FindAllMarkers(dat, only.pos = T)
markers_sig <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

write.csv(markers_sig, "3.Output/cd8/2.top20SigMarkerRes0.3.csv")
# Effector, Memory, and Dysfunctional CD8+ T Cell Fates in the Antitumor Immune Response
## https://www.hindawi.com/journals/jir/2016/8941260/tab1/

# effector 0
VlnPlot(dat, features = c("Klrg1", 
                          "Spn", # Cd43 
                          "Cd69", 
                          "Fasl", # Cd95
                          "Tnfrsf9" #Cd137
                          # "Tbx21",
                          # "Eomes",
                          # "Prdm1",
                          # "Runx3",
                          # "Stat4",
                          # "Stat5a",
                          # "Id2"
                          ),
        pt.size = 0,
        ncol = 5)

# central memory 4
VlnPlot(dat, features = c("Ccr7", 
                          "Cd44", # Cd43 
                          "Ptprc", # CD45
                          "Sell", # CD62L
                          "Il2rb", #CD122
                          "Il7r", # CD127
                          "Il2rb"#, # IL-15R
                          # "Tbx21",
                          # "Eomes",
                          # "Bcl6",
                          # "Tcf7",
                          # "Stat3",
                          # "Id3"
),
pt.size = 0,
ncol = 5)

# effector memory
VlnPlot(dat, features = c("Cd44",
                          "Ptprc",
                          "Il7r",
                          "Klrg1",
                          "Ccr7",
                          "Sell"
  
),
pt.size = 0,
ncol = 5)

# exhausted 2, 3
VlnPlot(dat, features = c("Ptprc",
                          "B3gat1",
                          "Fasl",
                          "Pdcd1",
                          "Ctla4",
                          "Havcr2",
                          "Lag3",
                          "Btla"
  
),
pt.size = 0,
ncol = 5)

# anergic/tolerant
VlnPlot(dat, features = c("Pdcd1",
                          "Lag3"
  
),
pt.size = 0,
ncol = 5)

# Senescent/regulatory
VlnPlot(dat, features = c("Klrg1",
                          "B3gat1",
                          "Cd28",
                          "Foxp3"
  
),
pt.size = 0,
ncol = 5)
