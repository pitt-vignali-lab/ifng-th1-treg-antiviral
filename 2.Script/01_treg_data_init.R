# ----------------------------------------------------------------------------
# Setup | Treg scRNA-seq QC, normalisation and initial clustering
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/1.datInit.R
# Input    : 1.Data/treg/ (CellRanger matrices, GSE223210)
# Output   : 1.Data/treg/TregDat.RDS
#            3.Output/treg/QC_*.png
# Note     : CellRanger 6.0.0, mm10
#            hashtag demultiplexing of individual mice
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(clustree)

# gene expression matrix
rawDat <- Read10X(paste0("1.Data/treg/filtered_feature_bc_matrix/"))
colnames(rawDat) <- gsub("-1", "", colnames(rawDat))

# hashtag matix
rawCite <- Read10X(paste0("1.Data/treg/umi_count/"), gene.column =1)
rawCite <- rawCite[grep("m_", rownames(rawCite))[1:8],]

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
cellsUse <- intersect(colnames(rawDat), colnames(rawCite))

# Subset RNA and HTO counts by joint cell barcodes
rawDat <- rawDat[, cellsUse]
rawCite <- as.matrix(rawCite[, cellsUse])

# Confirm that the HTO have the correct names
rownames(rawCite) <- paste0("HTO", "_", 1:8)

# Setup Seurat object
tagDat <- CreateSeuratObject(counts = rawDat)

# Normalize RNA data with log normalization
tagDat <- NormalizeData(tagDat)

# Find and scale variable features
tagDat <- FindVariableFeatures(tagDat, selection.method = "mean.var.plot")
tagDat <- ScaleData(tagDat, features = VariableFeatures(tagDat))

# Add HTO data as a new assay independent from RNA
tagDat[["HTO"]] <- CreateAssayObject(counts = rawCite)

# Normalize HTO data with centered log-ratio (CLR) transformation
tagDat <- NormalizeData(tagDat, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid clustering function for large applications You can also play with additional parameters (see documentation for HTODemux()) to adjust the threshold for classification Here we are using the default settings
tagDat <- HTODemux(tagDat, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(tagDat$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(tagDat) <- "HTO_maxID"
RidgePlot(tagDat, assay = "HTO", features = rownames(tagDat[["HTO"]])[1:3], ncol = 3)

# Compare number of UMIs for singlets, doublets and negative cells
Idents(tagDat) <- "HTO_classification.global"
VlnPlot(tagDat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Remove negative cells from the object
tagSub <- subset(tagDat, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(tagSub) <- "HTO"
tagSub <- ScaleData(tagSub, features = rownames(tagSub), 
                    verbose = FALSE)
tagSub <- RunPCA(tagSub, features = rownames(tagSub), approx = FALSE)
tagSub <- RunTSNE(tagSub, dims = 1:8, perplexity = 100)
Idents(tagSub) <-'HTO_classification'
Idents(tagSub) <-'HTO_classification.global'
DimPlot(tagSub)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(tagDat, assay = "HTO", ncells = 5000)

# Extract the singlets
allDat <- subset(tagDat, idents = "Singlet")

# Add tcr data ----
tcr <- read.csv("1.Data/treg/TCR/filtered_contig_annotations.csv")
tcr$barcode <- gsub("-1", "", tcr$barcode)
tcr <- tcr[!duplicated(tcr$barcode), ]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
rownames(tcr) <- tcr$barcode
tcr <- tcr[,-1]

allDat <- AddMetaData(object = allDat, metadata =tcr)
# ----

metaTreat <- data.frame(hash.ID = sort(as.character(unique(allDat@meta.data$hash.ID))),
                        reps  = c("UC1", paste0(rep("C",3),1:3),
                                       "UG1", paste0(rep("G",3),1:3)),
                        samples = c("UC", rep("C", 3),
                                         "UG", rep("G", 3)))
metaTreat <- allDat@meta.data["hash.ID"] %>%
  rownames_to_column("cells") %>%
  left_join(metaTreat, by = "hash.ID") %>%
  column_to_rownames("cells") %>%
  select(-1)


allDat <- AddMetaData(allDat, metaTreat)

# allDat$reps <- allDat$tubeLabel

rm(list = ls()[which(ls() != "allDat")])

nrow(allDat@meta.data)

table(allDat$reps)

DefaultAssay(allDat) <- "RNA"

allDat$log10GenesPerUMI <- log10(allDat$nFeature_RNA) / log10(allDat$nCount_RNA)
allDat$mitoRatio <- PercentageFeatureSet(object = allDat, pattern = "^mt-")/100

# Visualize the number of cell counts per sample ----
# png("3.Output/treg/QC_1__cell_number_distr.png", width = 8, height = 5, units = "in", res = 300)
allDat@meta.data %>% 
  ggplot(aes(x=reps, fill=reps)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# dev.off()

# Visualize the number UMIs/transcripts per cell ----
#  png("3.Output/treg/QC_2_cell_UMIs_distr.png", width = 8, height = 5, units = "in", res = 300)
allDat@meta.data %>% 
  ggplot(aes(color=reps, x=nCount_RNA, fill= reps)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA")
# dev.off()

# Visualize the distribution of genes detected per cell ----
# png("3.Output/treg/QC_3_genes_per_cell_distr.png", width = 8, height = 5, units = "in", res = 300)
allDat@meta.data %>% 
  ggplot(aes(color=reps, x=nFeature_RNA, fill= reps)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")
# dev.off()

allDat@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= reps, fill= reps)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs ----
allDat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~reps) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA vs. nFeature_RNA")

# Visualize the distribution of mitochondrial gene expression detected per cell ----
# png("3.Output/treg/QC_4_mitochondrial_gene_per_cell.png", width = 8, height = 5, units = "in", res = 300)
allDat@meta.data %>% 
  ggplot(aes(color=reps, x=mitoRatio, fill=reps)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.5) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Mitochondrial counts ratio")
# dev.off()

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI ----
allDat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = reps, fill=reps)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.81) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity")

# Filter out low quality reads using selected thresholds - these will change with experiment ----
filteredDat <- subset(x = allDat, 
                      subset= nCount_RNA >= 300 & 
                        nFeature_RNA >= 300 &
                        log10GenesPerUMI > 0.80 &
                        mitoRatio < 0.2)

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filteredDat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keepGenes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filteredCts <- counts[keepGenes, ]

# Reassign to filtered Seurat object
filteredDat <- filteredDat[, colnames(filteredDat) %in% colnames(filteredCts)]

rm(counts, nonzero, keepGenes, filteredCts)

filteredDat <- FindVariableFeatures(filteredDat, selection.method = "mean.var.plot")

# allGenes <- rownames(filteredDat)
filteredDat <- ScaleData(filteredDat, features = VariableFeatures(filteredDat))
# rm(allGenes)

filteredDat <- RunPCA(filteredDat, features = VariableFeatures(filteredDat))

# png("3.Output/treg/Analysis_1_PCA_in_every_sample.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(filteredDat, reduction = "pca", split.by = "reps")
# dev.off()

# Plot the elbow plot
# png("3.Output/treg/Analysis_2_ElbowPlot.png", width = 15, height = 5, units = "in", res = 300)
ElbowPlot(filteredDat, ndims = 50)
# dev.off()

# Determine percent of variation associated with each PC
pct <- filteredDat[["pca"]]@stdev / sum(filteredDat[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
# png("3.Output/treg/Analysis_3_ElbowPlot_PC_selection.png", width = 15, height = 5, units = "in", res = 300)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
# dev.off()

rm(plot_df, co1, co2, cumu, pct)

filteredDat <- FindNeighbors(filteredDat, dims = 1:pcs)
filteredDat <- FindClusters(filteredDat, resolution = seq(0.1, 1, 0.1))

filteredDat <- RunUMAP(filteredDat, dims = 1:pcs)
# filteredDat <- RunTSNE(filteredDat, dims = 1:pcs)

clustree(filteredDat, prefix = "RNA_snn_res.")

Idents(filteredDat) <- "RNA_snn_res.0.4"
# png("3.Output/treg/Analysis_4.1_overall_UMAP.png", width = 5, height = 5, units = "in", res = 300)
DimPlot(filteredDat, reduction = "umap", label = T)
# dev.off()

DimPlot(filteredDat, 
        reduction = "tsne",
        cells.highlight = WhichCells(filteredDat, idents = 6))

# png("3.Output/treg/Analysis_4.2_UMAP_check_batch_effect.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(filteredDat, reduction = "umap", split.by = "reps")
# dev.off()

# png("3.Output/treg/Analysis_5.1_overall_TSNE.png", width = 5, height = 5, units = "in", res = 300)
DimPlot(filteredDat, reduction = "tsne")
# dev.off()

# png("3.Output/treg/Analysis_5.2_TSNE_check_batch_effect.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(filteredDat, reduction = "tsne", split.by = "reps")
# dev.off()

# png("3.Output/treg/Analysis_6.1_Foxp3_expression_FeaturePlot.png", width = 25, height = 5, units = "in", res = 300)
DefaultAssay(filteredDat) <- "RNA"
FeaturePlot(filteredDat, features = "Foxp3", reduction = "umap", split.by = "reps", label = T)
# dev.off()

# png("3.Output/treg/Analysis_6.2_Foxp3_expression_DotPlot.png", width = 5, height = 5, units = "in", res = 300)
DotPlot(filteredDat, features = "Foxp3")
# dev.off()


saveRDS(filteredDat, "1.Data/treg/TregDat.RDS")
