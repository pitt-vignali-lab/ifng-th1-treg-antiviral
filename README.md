# IFNγ-driven TH1-like Tregs during viral infection

Analysis code for:

> **IFNγ-induction of TH1-like regulatory T cells controls antiviral responses**
> *Nature Immunology* (2023).
> doi:[10.1038/s41590-023-01453-w](https://doi.org/10.1038/s41590-023-01453-w)
> · PMC[10224582](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10224582/)

This repository contains the single-cell RNA-seq analysis code behind the
computational panels of the paper. Flow cytometry, mouse infection, adoptive
transfer and suppression assays are described in the paper and have no code
here.

---

## What the analysis does

Tregs and tetramer<sup>+</sup> CD8<sup>+</sup> T cells were sorted from
`Foxp3Cre-YFP` control and `Ifngr1`<sup>L/L</sup>`Foxp3Cre-YFP` mice, uninfected
or infected with LCMV clone 13 for 16 days, and profiled by 10x single-cell
RNA-seq with hashtag multiplexing and V(D)J libraries.

- Treg clusters are compared between genotypes to show that loss of IFNγ
  receptor signalling blocks TH1-like polarisation and shifts Tregs towards a
  TH2-like profile.
- Tetramer<sup>+</sup> CD8<sup>+</sup> T cells from the same mice are analysed
  to show how TH1-like Tregs restrain CD8 effector and memory responses.
- A public LCMV dataset (GSE186898) is mined for cytokine–receptor pairs across
  infected immune cell types.

---

## Repository layout

```
1.Data/
  treg/       Treg dataset and derived objects
  cd8/        tetramer+ CD8 dataset and derived objects
  revision/   integrated object, public data, curated reference tables
2.Script/     analysis scripts, numbered in run order
3.Output/     figures and tables produced by the scripts
```

Scripts expect the **repository root** as the R working directory:

```r
setwd("/path/to/ifng-th1-treg-antiviral")
source("2.Script/01_treg_data_init.R")
```

Each script carries a header giving its figure panel, inputs, outputs, the
upstream script it depends on, and its original filename in the lab archive.

The scripts fall into three groups, which were run in this order:

1. **`01`–`18` — Treg analysis.** The original Treg dataset: QC, clustering,
   differential expression and enrichment, plus the scripts that render the
   published Treg panels.
2. **`19`–`29` — CD8 analysis.** The tetramer<sup>+</sup> CD8 dataset, analysed
   the same way.
3. **`30`–`43` — Revision analyses.** Joint batch-corrected re-analysis of both
   datasets with `simspec`, pseudobulk comparisons, and the GSE186898
   cytokine–receptor data mining.

---

## Figure → script map

### Treg analysis

| Figure | Analysis | Script |
| --- | --- | --- |
| — | Treg QC, normalisation, initial clustering | [`01_treg_data_init.R`](2.Script/01_treg_data_init.R) |
| — | Cluster identity assignment from markers | [`02_treg_cluster_labeling.R`](2.Script/02_treg_cluster_labeling.R) |
| — | UMAP embedding and cluster overview | [`03_treg_umap.R`](2.Script/03_treg_umap.R) |
| — | Differential expression, control vs *Ifngr1*-deficient | [`04_treg_de_mast.R`](2.Script/04_treg_de_mast.R) |
| — | Differential expression within Treg subtypes | [`05_treg_de_mast_subtypes.R`](2.Script/05_treg_de_mast_subtypes.R) |
| — | GSEA and SPIA pathway analysis | [`06_treg_gsea.R`](2.Script/06_treg_gsea.R) |
| — | GSEA enrichment plots | [`07_treg_gsea_plots.R`](2.Script/07_treg_gsea_plots.R) |
| **3** | TH1 / TH2 marker feature and violin plots | [`08_treg_th1_th2_markers.R`](2.Script/08_treg_th1_th2_markers.R) |
| **2a** | UMAP of Treg clusters | [`09_fig2a_treg_umap.R`](2.Script/09_fig2a_treg_umap.R) |
| **2b** | Cluster distribution by genotype | [`10_fig2b_cluster_distribution.R`](2.Script/10_fig2b_cluster_distribution.R) |
| **2c** | Volcano plot of differentially expressed genes | [`11_fig2c_volcano.R`](2.Script/11_fig2c_volcano.R) |
| **3a** | Heat map of TH1 / TH2 gene expression | [`12_fig3a_th1_th2_heatmap.R`](2.Script/12_fig3a_th1_th2_heatmap.R) |
| **3b** | GSEA of immunologic signature (C7) gene sets | [`13_fig3b_gsea_c7.R`](2.Script/13_fig3b_gsea_c7.R) |
| **5h** | Heat map of cytokine, receptor, migratory and suppressive genes | [`14_fig5h_cytokine_receptor_heatmap.R`](2.Script/14_fig5h_cytokine_receptor_heatmap.R) |
| **ED 3e** | Cluster distribution | [`15_edfig3e_cluster_distribution.R`](2.Script/15_edfig3e_cluster_distribution.R) |
| **ED 3f** | Heat map of significant DE genes | [`16_edfig3f_de_heatmap.R`](2.Script/16_edfig3f_de_heatmap.R) |
| **ED 3g** | GSEA of the genotype comparison | [`17_edfig3g_gsea.R`](2.Script/17_edfig3g_gsea.R) |
| **ED 5a** | Foxp3-induced and Foxp3-repressed gene sets | [`18_edfig5a_foxp3_gene_sets.R`](2.Script/18_edfig5a_foxp3_gene_sets.R) |

### Tetramer<sup>+</sup> CD8 analysis

| Figure | Analysis | Script |
| --- | --- | --- |
| — | CD8 QC, normalisation, initial clustering | [`19_cd8_data_init.R`](2.Script/19_cd8_data_init.R) |
| **ED 9a** | Cluster identity assignment from markers | [`20_cd8_cluster_labeling.R`](2.Script/20_cd8_cluster_labeling.R) |
| **ED 9b–c** | Cluster composition statistics | [`21_cd8_cluster_statistics.R`](2.Script/21_cd8_cluster_statistics.R) |
| **ED 9a** | UMAP embedding | [`22_cd8_umap.R`](2.Script/22_cd8_umap.R) |
| — | Differential expression | [`23_cd8_de_analysis.R`](2.Script/23_cd8_de_analysis.R) |
| — | GSEA and SPIA pathway analysis | [`24_cd8_gsea.R`](2.Script/24_cd8_gsea.R) |
| — | GSEA enrichment and ridge plots | [`25_cd8_gsea_plots.R`](2.Script/25_cd8_gsea_plots.R) |
| **ED 9a** | UMAP of CD8 clusters | [`26_edfig9a_cd8_umap.R`](2.Script/26_edfig9a_cd8_umap.R) |
| **ED 9b–c** | Cluster distribution by genotype | [`27_edfig9bc_cluster_distribution.R`](2.Script/27_edfig9bc_cluster_distribution.R) |
| **ED 9f–g, 9k–n** | *Ifngr1* and migratory gene expression | [`28_edfig9_gene_expression.R`](2.Script/28_edfig9_gene_expression.R) |
| **ED 9** | Galaxy density plot of clusters | [`29_cd8_galaxy_plot.R`](2.Script/29_cd8_galaxy_plot.R) |

### Revision analyses

| Analysis | Script |
| --- | --- |
| Load Treg and CD8 datasets for joint analysis | [`30_rev_data_init.R`](2.Script/30_rev_data_init.R) |
| Batch correction with `simspec` and re-clustering | [`31_rev_batch_correction_simspec.R`](2.Script/31_rev_batch_correction_simspec.R) |
| Selected gene expression on the integrated UMAP | [`32_rev_selected_genes_umap.R`](2.Script/32_rev_selected_genes_umap.R) |
| Integrated UMAP and galaxy plot | [`33_rev_umap_galaxy.R`](2.Script/33_rev_umap_galaxy.R) |
| Differential abundance of clusters | [`34_rev_cluster_proportions.R`](2.Script/34_rev_cluster_proportions.R) |
| Pseudobulk aggregation | [`35_rev_pseudobulk.R`](2.Script/35_rev_pseudobulk.R) |
| Pseudobulk heat maps and volcano plots | [`36_rev_heatmaps_volcano.R`](2.Script/36_rev_heatmaps_volcano.R) |
| GSEA and SPIA on the pseudobulk comparison | [`37_rev_gsea.R`](2.Script/37_rev_gsea.R) |
| GSEA enrichment plots | [`38_rev_gsea_plots.R`](2.Script/38_rev_gsea_plots.R) |
| Load and preprocess public LCMV data (GSE186898) | [`39_gse186898_data_init.R`](2.Script/39_gse186898_data_init.R) |
| Annotate immune cell types in GSE186898 | [`40_gse186898_cluster_annotation.R`](2.Script/40_gse186898_cluster_annotation.R) |
| Cytokine–receptor pair interactions (`celltalker`) | [`41_gse186898_celltalker.R`](2.Script/41_gse186898_celltalker.R) |
| Circos plot of significant pairs | [`42_gse186898_celltalker_circos.R`](2.Script/42_gse186898_celltalker_circos.R) |
| Linkage view of the interaction network | [`43_gse186898_celltalker_linkage.R`](2.Script/43_gse186898_celltalker_linkage.R) |

---

## Data availability

**No sequencing data or intermediate objects are included in this repository.**

| Dataset | Source |
| --- | --- |
| Treg and tetramer<sup>+</sup> CD8 scRNA-seq generated for this study | GEO [**GSE223210**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223210) |
| Public LCMV scRNA-seq used for cytokine–receptor mining | GEO [GSE186898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186898) |

Download the CellRanger output for GSE223210 into `1.Data/treg/` and
`1.Data/cd8/`, and the GSE186898 matrix into
`1.Data/revision/GSE186898Data/`.

### Reference tables included here

`1.Data/revision/` carries the small curated tables that the cytokine–receptor
analysis depends on. They encode a manual curation step that cannot be
regenerated automatically, so they ship with the code:

| File | Contents |
| --- | --- |
| `mouseLRModified.csv` | Mouse ligand–receptor pairs, curated |
| `mouseLRFromPMID17201934Modified.csv` | Ligand–receptor pairs from PMID 17201934 |
| `mouseChemokine.csv`, `mouseChemokineReceptor.csv`, `mouseChemokineReceptorModified.csv` | Mouse chemokines and their receptors |
| `chemokinListFromPMID17201934.csv` | Chemokine list from PMID 17201934 |

`37_rev_gsea.R` additionally needs a local KEGG data directory for `SPIA` at
`1.Data/revision/RefData/`, generated with `SPIA::makeSPIAdata()`.

---

## Software

| Tool | Version / notes |
| --- | --- |
| CellRanger | 6.0.0, reference `refdata-gex-mm10-2020-A` |
| Seurat | standard workflow |
| simspec | `cluster_sim_spectrum`, `spectrum_type = corr_ztransform`, `corr_method = spearman` |
| presto | `wilcoxauc` for differential expression |
| MAST | via `Seurat::FindMarkers(test.use = "MAST")` |
| fgsea / msigdbr | MSigDB C2 and C7, set size 15–500 |
| clusterProfiler | customised `gseaplot` / `gseaplot2` |
| SPIA | KEGG pathway topology analysis |
| ComplexHeatmap | heat maps, batches scaled separately |
| celltalker | cytokine–receptor pair inference |
| ggraph | circos / network layouts |

Genes with adjusted *P* < 0.05 were taken as significant throughout.

---

## Contact

**Jian Cui** — Senior Bioinformatics Scientist  
Vignali Lab, Department of Immunology  
School of Medicine, University of Pittsburgh  
<cuij@pitt.edu>

Questions about the analysis or this code are welcome.

## License

[MIT](LICENSE)
