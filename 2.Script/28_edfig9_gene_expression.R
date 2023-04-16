# ----------------------------------------------------------------------------
# ED Fig 9f-g, 9k-n | Expression of Ifngr1 and migratory genes in CD8 T cells
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/7C.R
# Input    : 1.Data/cd8/cd8Dat.RDS
# Output   : 3.Output/cd8/7C.png
# Upstream : 20_cd8_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
require(scales)
library(ggpubr)
library(extrafont)
library(ggprism)

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


data <- data.frame(prop.table(table(Idents(datLab), datLab$reps),
                              margin = 2))
data$ifngr1 <- "indeficient"
data$ifngr1[grep("G", data$Var2)] <- "deficient"

# check the expression ----
AverageExpression(datLab, features = "Ifngr1", group.by = "reps", assays = "RNA")[[1]]

AverageExpression(subset(datLab, idents = "Memory"), features = "Ifngr1", group.by = "reps", assays = "RNA")[[1]]


## Bar chart
dataBar <- data[data$Var1 == "Memory",] %>% 
  group_by(ifngr1) %>% 
  summarise(mean = mean(Freq),
            sd = sd(Freq))

png("3.Output/cd8/7C.png", width = 3, height = 8, units = "in", res = 300)
ggplot(data[data$Var1 == "Memory",], aes(fill = ifngr1))  +
  # stat_compare_means(aes(x = ifngr1, y = Freq*100, group = ifngr1),
  #                    method = "anova",
  #                    label =  "p",
  #                    label.x = 1.5,
  #                    label.y = 5.5, 
  #                    size = 10) + 
  stat_summary(data = dataBar, aes(x = ifngr1, y = mean*100), na.rm = TRUE,
               geom = "col", fun = mean, colour = "black", size = 1.5, alpha =0.7) +
  geom_point(aes(x = ifngr1, y = Freq*100),col = "black",shape = 21,
             position=position_jitterdodge(jitter.width = 0.8),
             size = 12, stroke = 2)+
  scale_fill_manual(name = "Status",
                     breaks = c("deficient", "indeficient"),
                     values=c("#FF0000", "#7434a4"),
                     label = c(expression(paste(UP," ","in"," ",Foxp3^{Cre-YFP})),
                               expression(paste(UP," ","in"," ",Ifngr1^{L/L},Foxp3^{Cre-YFP})))) +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 5.5), expand = c(0,0)) +
  ylab("% in memory cluster")+
  xlab(NULL)+
  theme_prism() +
  theme(text = element_text(family = "Arial", size = 15, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 40), 
        axis.title.y = element_text(size = 40),
        legend.position = "bottom")+
  # ggtitle("Effects of Ifngr1 mutation") +
  NoLegend()
dev.off()
