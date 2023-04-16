# ----------------------------------------------------------------------------
# Data mining | Load and preprocess the public LCMV dataset GSE186898
#
# Paper    : IFNgamma-induction of TH1-like regulatory T cells controls
#            antiviral responses
# Original : 2.Script/6.3.1.gse186898DatInit.R
# Input    : 1.Data/revision/GSE186898Data/ (download from GEO)
# Output   : 1.Data/revision/
# Note     : External dataset, not generated in this study
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/")){
  dir.create("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/")
}

library(tidyverse)
library(Seurat)

## gene expression matrix
dat <- Read10X(data.dir = "1.Data/revision/GSE186898Data/")

datRaw <- dat$`Gene Expression`

# hashtag
citeRaw <- dat$`Antibody Capture`[1:8,]

cellsUse <- intersect(colnames(datRaw), colnames(citeRaw))

# Subset RNA and HTO counts by joint cell barcodes
datRaw <- datRaw[, cellsUse]
citeRaw <- as.matrix(citeRaw[, cellsUse])

# Setup Seurat object
datTag <- CreateSeuratObject(counts = datRaw)

# Normalize RNA data with log normalization
datTag <- NormalizeData(datTag)

# Find and scale variable features
datTag <- FindVariableFeatures(datTag, selection.method = "mean.var.plot")
datTag <- ScaleData(datTag, features = VariableFeatures(datTag))

# Add HTO data as a new assay independent from RNA
datTag[["HTO"]] <- CreateAssayObject(counts = citeRaw)

# Normalize HTO data with centered log-ratio (CLR) transformation
datTag <- NormalizeData(datTag, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid clustering function for large applications You can also play with additional parameters (see documentation for HTODemux()) to adjust the threshold for classification Here we are using the default settings
datTag <- HTODemux(datTag, assay = "HTO", positive.quantile = 1-10E-8)

# Global classification results
table(datTag$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(datTag) <- "HTO_maxID"
RidgePlot(datTag, assay = "HTO", features = rownames(datTag[["HTO"]])[1:8], ncol = 4)+
  coord_fixed()
ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/ridgePlotofHTOs.png"), 
       width = 12, height = 8, units = "in", bg = "white")

# Compare number of UMIs for singlets, doublets and negative cells
Idents(datTag) <- "HTO_classification.global"
datTag$HTO_classification.global <- factor(datTag$HTO_classification.global,
                                           levels = c("Singlet", "Doublet", "Negative"))
VlnPlot(datTag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)+
  coord_fixed()
ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/barPlotofHTOs.png"), 
       width = 5, height = 4, bg = "white")

# Extract the singlets
dat <- subset(datTag, idents = "Singlet")

sampDf <- read.table("1.Data/revision/GSE186898Data/GSE186898_readme.txt", 
                     sep = "\t")[-c(1:2),]
colnames(sampDf) <- sampDf[1,]
sampDf <- sampDf[-1,]

sampDf$`Sample name` <- paste0("IL27-sample", 1:8)
names(sampDf)[1] <- "hash.ID"

dat <- AddMetaData(dat, 
                   metadata = dat@meta.data %>% select("hash.ID") %>% 
                     rownames_to_column(var = "cellID") %>% 
                     merge(sampDf, 
                           by = "hash.ID") %>% 
                     column_to_rownames(var = "cellID"))

dat$log10GenesPerUMI <- log10(dat$nFeature_RNA) / log10(dat$nCount_RNA)
dat$mitoRatio <- PercentageFeatureSet(object = dat, pattern = "^mt-")/100
dat$ribonRatio <- PercentageFeatureSet(dat, pattern = "^Rp[sl]")/100
dat$hbRatio <- PercentageFeatureSet(dat, pattern = "^Hb[^(p)]")/100

VlnPlot(dat, features = c("nFeature_RNA")) +
  scale_y_continuous(trans='log10') 

VlnPlot(dat, features = c("nCount_RNA")) +
  scale_y_continuous(trans='log10')

# cell number distribution
dat@meta.data %>% 
  ggplot(aes(x=title, fill=title)) + 
  geom_bar() +
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("NCells") +
  NoLegend()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/cellNumberDist.png"), 
       width = 5, height = 4, bg = "white")

# number of UMIs/transcripts per replicate
dat@meta.data %>% 
  ggplot(aes(color=title, x=nCount_RNA, fill= title)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA")

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/umiDist.png"), 
       width = 5, height = 4, bg = "white")

# genes detected per replicate
dat@meta.data %>% 
  ggplot(aes(color=title, x=nFeature_RNA, fill= title)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/geneDist.png"), 
       width = 5, height = 4, bg = "white")

dat@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= title, fill= title)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_hline(yintercept = log10(250)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("nFeature_RNA") +
  NoLegend()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/geneDistBox.png"), 
       width = 5, height = 4, bg = "white")

# the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
dat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 250) +
  facet_wrap(~title) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        text = element_text(size = 6)) +
  ggtitle("nCount_RNA vs. nFeature_RNA")

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/umiVsGenes.png"), 
       width = 5, height = 5, bg = "white")

# mitochondrial gene expression
dat@meta.data %>% 
  ggplot(aes(color=title, x=mitoRatio, fill=title)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial counts ratio")

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/mitochondrialGeneDist.png"), 
       width = 5, height = 4, bg = "white")

# verall complexity of the gene expression by visualizing the genes detected per UMI
dat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = title, fill=title)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity")

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/complexity.png"), 
       width = 5, height = 4, bg = "white")

# Filter out low quality reads using selected thresholds - these will change with experiment ----
filteredDat <- subset(x = dat, 
                      subset= nCount_RNA >= 300 & 
                        nFeature_RNA >= 250 &
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
filteredDat <- filteredDat[, colnames(filteredDat) %in% 
                             colnames(filteredCts)]

# cell number distribution
filteredDat@meta.data %>% 
  ggplot(aes(x=title, fill=title)) + 
  geom_bar() +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("NCells") +
  NoLegend()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/cellNumberDistPostFilt.png"), 
       width = 5, height = 4, bg = "white")

# dat processing
rm(counts, nonzero, keepGenes, filteredCts)

filteredDat <- FindVariableFeatures(filteredDat, 
                                    selection.method = "mean.var.plot")

# allGenes <- rownames(filteredDat)
filteredDat <- ScaleData(filteredDat, 
                         features = VariableFeatures(filteredDat))
# rm(allGenes)

filteredDat <- RunPCA(filteredDat, 
                      features = VariableFeatures(filteredDat))

DimPlot(filteredDat, reduction = "pca") +
  coord_fixed() +
  NoLegend()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/pcaOverall.png"), 
       width = 4, height = 5, bg = "white")

DimPlot(filteredDat, reduction = "pca", split.by = "title") +
  coord_fixed() +
  NoLegend()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/pcaRep.png"), 
       width = 16, height = 4, bg = "white")

ElbowPlot(filteredDat, ndims = 50)

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/elbow.png"), 
       width = 15, height = 4, bg = "white")

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
# png("3.Output/revision/Analysis_3_ElbowPlot_PC_selection.png", width = 15, height = 5, units = "in", res = 300)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

ggsave(paste0("3.Output/revision/6.Chemokins/GSE186898/1.Initialization/elbowPCselect.png"), 
       width = 15, height = 4, bg = "white")

rm(plot_df, co1, co2, cumu, pct)

saveRDS(filteredDat, "3.Output/revision/6.Chemokins/GSE186898/1.initDat.RDS")
