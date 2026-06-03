# PTM-Placenta-Atlas

Link to the bioRxiv manuscript: https://www.biorxiv.org/content/10.1101/2025.08.15.669966v1

Link to the GEO repository of FASTQ files: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305531

## Overview
This repository contains the code and instructions to reproduce the results from:
> **"A single-cell transcriptomic atlas of the pigtail macaque placenta in late gestation "**, Amanda Li, Richard Li, Hazel Huang, Hong Zhao, Briana Del Rosario, Miranda Li, Edmunda Li, Andrew Vo, Gygeria Manuel, Orlando Cervantes, Raj Kapur, Jeff Munson, Austyn Orvis, Michelle Coleman, Melissa Berg, Britni Curtis, Brenna Menz, Jin Dai, Inah Golez, Solomon Wangari, Chris English, Audrey Baldessari, Lakshmi Rajagopal, John Cornelius, Kristina Adams Waldorf 

## Repository Structure
- `src/` — Source code.


# Macaque Placenta Single-Cell Atlas: Analysis Pipeline

---

## Scripts

| Script | Language | Purpose |
|---|---|---|
| `placenta_v3_organized.R` | R | QC, integration, annotation, compositional analysis, figures |
| `run_scvelo_full.py` | Python | RNA velocity: full atlas (stochastic model) |
| `run_scvelo_v2.py` | Python | RNA velocity: cell type subsets (dynamical model, stochastic fallback) |

---

## Pipeline Steps

**`placenta_v3_organized.R`**

1. **Cell calling & ambient RNA removal**: emptyDrops (FDR < 0.01; fallback to CellRanger filtered matrix if <100 cells), SoupX ambient correction per sample
2. **Doublet removal & QC**: scDblFinder; retain singlets with 200–7,500 genes and <2.5% MT reads
3. **Normalization & integration**: SCTransform (regress MT%), PCA (50 PCs), Harmony batch correction (theta = 2, 20 iterations, 30 PCs), UMAP (30 neighbors, min.dist = 0.3)
4. **Label transfer**: FindTransferAnchors / TransferData from reference atlas (SCT, 30 PCs)
5. **Subset annotation**: Lineages (macrophage, T/NK/B, trophoblast, stromal/endothelial) re-clustered independently (SCTransform → PCA → Harmony → UMAP, res = 0.15); clusters annotated by marker dot plots, interactive lasso gating, KNN-based reassignment, and FindSubCluster; mast cell contaminants (MS4A2/CPA3/FCER1A) removed; final labels propagated to full atlas (`cell_type_v5`, 38 cell types)
6. **Compositional analysis**: scComp Bayesian model comparing treatment (Media vs. Saline), tissue compartment, and fetal sex (FDR < 5%)
7. **Figure generation**: Atlas and subset UMAPs, lineage dot plots, tissue/sex composition bar plots

**`run_scvelo_full.py` and `run_scvelo_v2.py`**

1. **Load spliced/unspliced counts**: velocyto loom files concatenated; barcodes converted from `sampleName:BARCODEx` → `BARCODE-1` format
2. **Barcode matching**: Cells matched to Seurat metadata CSV via composite `sampleID_barcode` key; UMAP coordinates and annotations transferred to AnnData
3. **Preprocessing**: `scv.pp.filter_and_normalize` (min_shared_counts = 20, n_top_genes = 2,000); `scv.pp.moments` (n_pcs = 30, n_neighbors = 30)
4. **Velocity**: Full atlas uses stochastic model. Subsets (`run_scvelo_v2.py`) use dynamical model (`recover_dynamics`) with stochastic fallback if >90% NaN velocity values
5. **Outputs**: Stream plots (SVG), arrow plots (PDF), confidence/pseudotime/latent time scatter plots; top velocity genes per cell type (`rank_velocity_genes`, min_corr = 0.3)

---

## Inputs

| File | Description |
|---|---|
| `cell_ranger_outs/<sample>/outs/` | CellRanger output (raw + filtered matrices) |
| `placenta_seurat_v7.RDS` | Reference atlas for label transfer |
| `velocyto_out/*.loom` | Spliced/unspliced counts |
| `full_placenta_metadata.csv` | Seurat metadata + UMAP coords (full atlas) |

## Outputs

| File | Description |
|---|---|
| `placenta_seurat_FINAL_v4.RDS` | Final annotated Seurat object |
| `mac_subset_v5.RDS` / `tb_subset_FINAL.RDS` / `t_subset_v3.RDS` | Annotated lineage subsets |
| `full_dataset_velocity.h5ad` / `TB_velocity.h5ad` | scVelo results |
| `sccomp_*_FDR.csv` | Compositional analysis results |
| `top_20_markers_*.csv` | Top marker genes per cell type |

---

## Dependencies

**R:** `Seurat` (v5.1.0), `harmony` (v1.2.0), `DropletUtils`, `SoupX`, `scDblFinder`, `sccomp`, `scCustomize`, `SeuratExtend`, `dittoSeq`, `ggplot2`, `ggrepel`, `shadowtext`, `patchwork`, `plotly`, `shiny`, `FNN`, `tidyverse`, `future`

**Python:** `scvelo`, `scanpy`, `pandas`, `numpy`, `matplotlib`, `scipy`
