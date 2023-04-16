# ----------------------------------------------------------------------------
# ED Fig 9b-c | CD8 cluster distribution by genotype
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/7B.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/7B.png
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(extrafont)

loadfonts(device = "win", quiet = TRUE)

dat <- readRDS("1.Data/cd8/cd8Dat.RDS")
DefaultAssay(dat) <- "RNA"
Idents(dat) <- "RNA_snn_res.0.3"

datSub <- subset(dat, idents = 8, invert = TRUE)

clusterId <- c("Activated", # 0
               "Cytotoxic", # 1
               "Progenitor\n exhausted", # 2
               "Proliferating", # 3
               "Exhausted", # 4
               "Memory", # 5
               "IFN signaling", # 6
               "Metabolically\n hyperactive" # 7
)

names(clusterId) <- levels(datSub)

datLab <- RenameIdents(datSub, clusterId)

DimPlot(datLab, label = T) + NoLegend()

png("3.Output/cd8/7B.png", width = 9*(0.5^0.5), height = 4.5*(.5^0.5), units = "in", res = 300)
data.frame(prop.table(table(Idents(datLab), datLab$reps),
                      margin = 2)) %>% 
  mutate(Var1 = factor(paste0(Var1, ifelse(Var1 == "Memory", "**", "")),
                       levels = c("Activated", # 0
                                  "Cytotoxic", # 1
                                  "Progenitor\n exhausted", # 2
                                  "Proliferating", # 3
                                  "Exhausted", # 4
                                  "Memory**", # 5
                                  "IFN signaling", # 6
                                  "Metabolically\n hyperactive"))) %>%  
  mutate(Var3 = ifelse(Var2 %in% Var2[grep("U", Var2)], "Na", "Cl")) %>% 
  mutate(Var4 = ifelse(Var2 %in% Var2[grep("C", Var2)], "Fx", "Ifn")) %>% 
  ggplot(aes(x = Freq, y = Var2, 
             fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.7, 
           position = position_fill(reverse = TRUE)) +
  scale_x_continuous(labels = scales::percent) +
  ylab(NULL) +
  xlab("Percentage") +
  # facet_grid(Var3 + Var4 ~ ., scales = "free_y", space = "free_y") +
  theme_bw()+
  scale_fill_discrete(name = "Cell type") +
  theme(legend.position="bottom",
        text = element_text(family = "Arial", size = 15, face = "bold"))
dev.off()
