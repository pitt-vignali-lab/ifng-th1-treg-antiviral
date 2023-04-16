# ----------------------------------------------------------------------------
# ED Fig 3e | Cluster distribution of control versus Ifngr1-deficient Tregs
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/S2A.R
# Input    : 1.Data/treg/TregDat.RDS
# Output   : 3.Output/treg/S2A.png
# Upstream : 02_treg_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(ggh4x)
library(extrafont)

loadfonts(device = "win", quiet = TRUE)

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
               "Metabolically\n hyperactive", # 2
               "Proliferating", # 3
               "Gene transcription and\n protein stability" # 4
)

names(clusterId) <- levels(dat)

data.frame(table())

datLab <- RenameIdents(dat, clusterId)

DimPlot(datLab, label = T) + NoLegend()


png("3.Output/treg/S2A.png", width = 9*(0.75^0.5), height = 5*(0.75^0.5), units = "in", res = 300)
data.frame(prop.table(table(Idents(datLab), datLab$reps),
                      margin = 2)) %>% 
  mutate(Var1 = factor(paste0(Var1, ifelse(Var1 == "Metabolically\n hyperactive", "*", "")),
                       levels = c("Effector", # 0
                                  "Central", # 1
                                  "Metabolically\n hyperactive*", # 2
                                  "Proliferating", # 3
                                  "Gene transcription and\n protein stability"))) %>%  
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

