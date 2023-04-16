# ----------------------------------------------------------------------------
# Revision | Differential abundance of clusters between genotypes
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/3.1.DiffBehaClst(ProportionStatistic).R
# Input    : 1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS
# Output   : 3.Output/revision/3.DiffBehaClst/
# Upstream : 31_rev_batch_correction_simspec.R
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/3.DiffBehaClst")){
  dir.create("3.Output/revision/3.DiffBehaClst")
}

library(Seurat)
library(tidyverse)
library(ggpubr)

# method 1, cluster propertion statistical comparison ----
dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")
DimPlot(dat)

repVsClst <- data.frame(prop.table(table(dat$newAnno, dat$reps), margin = 1)) |> 
  rename("Cluster" = 1, "Replicate" = 2) |> 
  mutate(Sample = gsub("[0-9]", "", Replicate)) |> 
  mutate(Sample = factor(Sample, levels = c("UC", "UG", "C", "G"))) 

repVsClst$Infection <- "Infected"
repVsClst[grep("U", repVsClst$Sample),]$Infection <- "Uninfected"

repVsClst$Genotype <- "Indeficient"
repVsClst[grep("G", repVsClst$Sample),]$Genotype <- "Deficient"



# Allowed methods are one of: t.test, wilcox.test, anova, kruskal.test, kruskal.test
ggplot(repVsClst, aes(x = Infection, y = Freq)) + 
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~Cluster) + 
  labs(title = "Effects of LCMV infection") +
  stat_compare_means()

ggplot(repVsClst, aes(x = Genotype, y = Freq)) + 
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~Cluster) + 
  labs(title = "Effects of Ifngr1 deficiency")+
  stat_compare_means()

my_comparisons = list(c("UG", "UC"), c("C", "UC"), c("G", "UC"),
                      c("C", "UG"), c("G", "UG"), 
                      c("G", "C"))

ggplot(repVsClst, aes(x = Sample, y = Freq)) + 
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~Cluster) + 
  stat_compare_means(comparisons = my_comparisons)

# method 2.MiloR ----
library(miloR)
library(SingleCellExperiment)
library(scater)
library(patchwork)

dat <- readRDS("1.Data/revision/2.datIntegratedWithBatchEffectsRemoval.RDS")
shortAnno <- data.frame(newAnno = unique(dat$newAnno),
                        shortAnno = c("MetaHyp",
                                      "Eff",
                                      "Centr",
                                      "GeneTrans",
                                      "Prolf"))
dat <- AddMetaData(dat, merge(dat@meta.data["newAnno"] |> 
                                rownames_to_column(var = "cell"),
                              shortAnno, by = "newAnno") |> 
                     select("cell", "shortAnno") |> 
                     column_to_rownames(var = "cell"))



dat$idtp <- paste(dat$shortAnno, dat$samples, sep = "_")

datSce <- as.SingleCellExperiment(dat)

milo <- Milo(datSce)

milo

milo <- buildGraph(milo, k = 30, d = 20)

milo <- makeNhoods(milo, prop = 0.2, 
                   k = 30, 
                   d=20, 
                   refined = TRUE)

plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   sample="idtp")

head(nhoodCounts(milo))

# Defining experimental design 
designDf <- data.frame(colData(milo))[,c("shortAnno", "samples", "idtp")] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "idtp")

designDf

# Computing neighbourhood connectivity
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

# Testing
da_results <- testNhoods(milo, 
                         design = ~ shortAnno + samples, 
                         design.df = designDf)

da_results %>%
  arrange(SpatialFDR) %>%
  head()

# Inspecting DA testing results
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggsave("3.Output/revision/3.DiffBehaClst/uncorrectedPVal.png", 
       width = 5, height = 4, units = "in")

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05))
ggsave("3.Output/revision/3.DiffBehaClst/uncorrectedPVal(Volcano).png", 
       width = 5, height = 4, units = "in")

milo <- buildNhoodGraph(milo)

umap_pl <- plotUMAP(milo, 
                    colour_by="shortAnno", 
                    text_by = "shortAnno", 
                    text_size = 3, 
                    point_size=1) +
  guides(fill="none") +
  coord_fixed()

nh_graph_pl <- plotNhoodGraphDA(milo, 
                                da_results, 
                                alpha=0.5) +
  coord_fixed()


umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

ggsave("3.Output/revision/3.DiffBehaClst/DAplot.png", 
       width = 12, height = 7, unit = "in")


da_results <- annotateNhoods(milo, da_results, coldata_col = "newAnno")

head(da_results)

ggplot(da_results, aes(idtp_fraction)) + geom_histogram(bins=50)

da_results$celltype <- ifelse(da_results$idtp_fraction < 0.6, 
                              "Mixed", 
                              da_results$idtp)

plotDAbeeswarm(da_results, group.by = "idtp")
ggsave("3.Output/revision/3.DiffBehaClst/distDAChangeCommon.png",
       width = 8, height = 6, units = "in")

plotDAbeeswarm(da_results, group.by = "celltype")
ggsave("3.Output/revision/3.DiffBehaClst/distDAChange.png",
       width = 8, height = 6, units = "in")
# head(da_results)
# 
# ggplot(da_results, aes(ident_fraction)) + geom_histogram(bins=50)
# 
# da_results$celltype <- ifelse(da_results$ident_fraction < 0.7, 
#                               "Mixed", 
#                               da_results$ident)
# 
# plotDAbeeswarm(da_results, group.by = "celltype")

# Identifying signatures of DA subpopulations ----
logcounts(milo) <- log1p(counts(milo))

dge_smp <- findNhoodMarkers(milo, 
                            da_results,
                            assay = "counts", 
                            gene.offset = FALSE, 
                            da.fdr = 0.05,
                            aggregate.samples = TRUE, 
                            sample_col = "idtp",
                            subset.nhoods = da_results$ident %in% c(0))

head(dge_smp)

markers <- dge_smp[which(dge_smp$adj.P.Val_1 < 0.05 ), "GeneID"]
logcounts(milo) <- log1p(counts(milo))
milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionDA(milo, da_results, features = markers,
                      subset.nhoods = da_results$ident %in% c(0),
                      assay="logcounts",
                      scale_to_1 = TRUE, cluster_features = TRUE)






DotPlot(dat, features = markers, cols = c("white", "blue"), col.min = 0)

# Finding markers of DA populations ----

## Add log normalized count to Milo object
milo <- logNormCounts(milo)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(milo, da_results, 
                                          subset.row = rownames(milo)[1:10])

head(da_nhood_markers)

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo)[1:10], 
                                          aggregate.samples = TRUE, 
                                          sample_col = "idtp")

head(da_nhood_markers)

## Run buildNhoodGraph to store nhood adjacency matrix
milo <- buildNhoodGraph(milo)

## Find groups
da_results <- groupNhoods(milo, da_results, max.lfc.delta = 10)
head(da_results)

plotNhoodGroups(milo, da_results, layout="UMAP")

plotDAbeeswarm(da_results, "NhoodGroup")

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 5, overlap=5)
plotNhoodGroups(milo, da_results, layout="UMAP")

plotDAbeeswarm(da_results, group.by = "NhoodGroup")
