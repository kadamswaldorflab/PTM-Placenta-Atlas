setwd("/your/working/directory")
library(Seurat)
library(tidyverse)
library(data.table)
library(scCustomize)
library(sccomp)
library(speckle)
library(loupeR)
library(forcats)
library(ggplot2)
library(magrittr)
library(patchwork)
library(openxlsx)
library(dplyr)
library(ggrepel)
library(SeuratExtend)
library(monocle3)
library(SeuratDisk)
library(Seurat.utils)
library(scanalysis)
library(loomR)
library(SeuratWrappers)
library(convert2anndata)
library(anndata)

set.seed = 42
placenta_seurat <- readRDS("placenta_seurat.RDS")
immune_sub <- readRDS("immune_sub_v9.RDS")
macro_sub <- readRDS("macro_sub_v8.RDS")
placenta_sce <- readRDS("placenta_sce.RDS")

# Get cell counts per sample
meta_data <- placenta_seurat@meta.data %>% as.data.table
count_table <- meta_data[, .N, by = "sample_id"]
count_table
write.csv(count_table, file = "count_table.csv")

# Percent expressing for feature of interest
percent <- Percent_Expressing(seurat_object = placenta_seurat, features = c("NOTCH2"))
percent

percent <- t(percent) # convert columns to rows
percent <- as.data.frame(percent) # need to change it to a data frame as t() converts it to a matrix
percent <- rownames_to_column(percent, "cluster")
# percent[c('cluster', 'NOTCH2')] <- str_split_fixed(percent$cell_type, "_", 2)
# percent <- percent[c("cell_type", "p_stat", "segment_1")]
write.csv(percent, "notch2_cluster_cell_percentages.csv")
ggplot(data=percent, aes(x=cluster, y=NOTCH2)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  labs(title = "NOTCH2 Cell Percentage", x= "Cluster", y= "NOTCH2 Percent") + 
  guides(fill = guide_legend(reverse=TRUE))

  #### All markers ####

all_markers <- FindAllMarkers(object = placenta_seurat) %>%
  Add_Pct_Diff()

write_csv(all_markers, file = "all_markers_placenta.csv")

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")

plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers, flip = TRUE, x_lab_rotate = 90)
plots[[1]]

clusters <- all_markers[all_markers["cluster"]==c("6", "12", "17", "20"),]

top_markers <- Extract_Top_Markers(marker_dataframe = clusters, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")

placenta_seurat[["seurat_clusters"]] <- Idents(placenta_seurat);
clustersToSubset <- c("6", "12", "17", "20");
data.sub <- subset(placenta_seurat, idents = clustersToSubset)

plots <- Clustered_DotPlot(seurat_object = data.sub, features = top_markers, flip = TRUE, x_lab_rotate = 90)
plots[[1]]

# Get cell counts per cluster per tissue
cell_counts <- table(placenta_seurat@meta.data$seurat_clusters, placenta_seurat@meta.data$tissue_type)
cell_counts <- as.data.frame(cell_counts)
cell_counts <- rownames_to_column(cell_counts, "cluster")


ggplot(data=cell_counts, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill") + 
  labs(title = "Cell Counts Per Cluster Per Tissue", x= "Cluster", y= "Cell Counts") + 
  guides(fill = guide_legend(title = "Tissue",reverse=FALSE)) +
  scale_fill_discrete(labels = c("CV", "Decidua", "Meminoc"))

  # Immune cell subset
Idents(placenta_seurat) <- "seurat_clusters"

clustersToSubset <- c("1", "8", "29", "7", "26", "27", "31", "23")
data.sub <- subset(placenta_seurat, idents = clustersToSubset)

DimPlot_scCustom(data.sub)

# Run PCA and UMAP
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

DimPlot_scCustom(data.sub)
Idents(data.sub) <- "celltype18"

DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 8,
                                                                                                                                                          palette = "varibow", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 8, palette = "ditto_seq", shuffle_pal = FALSE), color_seed = 5)

saveRDS(data.sub, file = "immune_sub_v8.RDS")


# Macro subset

Idents(placenta_seurat) <- "seurat_clusters"

clustersToSubset <- c("9", "3", "15", "2", "14", "0", "10")
data.sub <- subset(placenta_seurat, idents = clustersToSubset)

DimPlot_scCustom(data.sub)

# Run PCA and UMAP
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.4)

DimPlot_scCustom(data.sub)
Idents(data.sub) <- "celltype18"

DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 8,
                                                                                                                                                   palette = "varibow", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 8, palette = "ditto_seq", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 8, palette = "polychrome", shuffle_pal = FALSE), color_seed = 5)

saveRDS(data.sub, file = "macro_sub_v8.RDS")

# Convert immune to anndata for scVelo analysis using SeuratExtend
Seu2Adata(immune_sub)
# Save the AnnData object to a local file
immune_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "immune_adata.h5ad")
adata.Save(immune_adata_path)

# Convert macro to anndata for scVelo analysis using SeuratExtend
Seu2Adata(macro_sub)
# Save the AnnData object to a local file
macro_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "macro_adata.h5ad")
adata.Save(macro_adata_path)

# CTB_EC_STB subset

Idents(placenta_seurat) <- "celltype19"
clustersToSubset <- c("STB", "PC", "EC1", "EC2", "FIB", "CTB", "iCTB", "schCTB")
data.sub <- subset(placenta_seurat, idents = clustersToSubset)

DimPlot_scCustom(data.sub)

# Run PCA and UMAP
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.4)

DimPlot_scCustom(data.sub)
Idents(data.sub) <- "celltype19"

DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9,
                                                                                                                                                   palette = "varibow", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9, palette = "ditto_seq", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9, palette = "polychrome", shuffle_pal = FALSE), color_seed = 5)

saveRDS(data.sub, file = "CTB_EC_STB_v2.RDS")
data.sub <- readRDS("CTB_EC_STB_v2.RDS")

# Convert to anndata for scVelo analysis using SeuratExtend
Seu2Adata(data.sub)
# Save the AnnData object to a local file
CTB_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "CTB_EC_STB_adata.h5ad")
adata.Save(CTB_adata_path)


# DSC subset

Idents(placenta_seurat) <- "celltype19"
clustersToSubset <- c("DMSC", "DSC0", "DSC1", "DSC2", "DSC3", "DSC4")
data.sub <- subset(placenta_seurat, idents = clustersToSubset)

DimPlot_scCustom(data.sub)

# Run PCA and UMAP
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.4)

DimPlot_scCustom(data.sub)
Idents(data.sub) <- "celltype19"

DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9,
                                                                                                                                                   palette = "varibow", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9, palette = "ditto_seq", shuffle_pal = FALSE), color_seed = 5)
DimPlot_scCustom(data.sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9, palette = "polychrome", shuffle_pal = FALSE), color_seed = 5)

saveRDS(data.sub, file = "DSC_v2.RDS")


# Convert to anndata for scVelo analysis using SeuratExtend
Seu2Adata(data.sub)
# Save the AnnData object to a local file
DSC_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "DSC_adata.h5ad")
adata.Save(DSC_adata_path)


# Subset IgM+B cluster to find MAST cells
Graphs(immune_sub)
immune_sub <- FindSubCluster(immune_sub, cluster = "IgM+B", graph.name = "SCT_snn", subcluster.name = "mast", resolution = 0.1)
DimPlot_scCustom(immune_sub, group.by = "mast", label = TRUE)

# renaming immune subclusters
Idents(immune_sub) <- "mast"
levels(immune_sub)

new.cluster.ids <- c("TRM CD8+T", "uNK1" ,     "uNK2" ,     "NKT" ,      "CD4+Treg" , "IgG+B",     "uNKP" ,     "IgM+B" ,  "Mast",  
                     "IgM+B")
names(new.cluster.ids) <- levels(immune_sub)
immune_sub <- RenameIdents(immune_sub, new.cluster.ids)
DimPlot_scCustom(immune_sub, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 10,
                                                                                                                                     palette = "polychrome", shuffle_pal = TRUE), color_seed = 42)
immune_sub$merged <- Idents(immune_sub)
saveRDS(immune_sub, file = "immune_sub_v9.RDS")

# Add immune cell labels back into main UMAP
CellsMetaTrim <- subset(immune_sub@meta.data, select = c("merged"))
placenta_seurat <- AddMetaData(placenta_seurat, CellsMetaTrim)

meta.data <- placenta_seurat@meta.data

meta.data <- meta.data %>% mutate(celltype19 = coalesce(merged, celltype18))
celltype19 <- subset(meta.data, select = c("celltype19"))
placenta_seurat <- AddMetaData(placenta_seurat, celltype19)

placenta_seurat$merged <- NULL

Idents(placenta_seurat) <- "celltype19"

saveRDS(placenta_seurat, file = "placenta_seurat_v20.RDS")

# Convert placenta seurat to anndata for scVelo analysis using SeuratExtend
Seu2Adata(placenta_seurat)
# Save the AnnData object to a local file
placenta_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "placenta_adata.h5ad")
adata.Save(placenta_adata_path)

# Convert immune seurat to anndata for scVelo analysis using SeuratExtend
Seu2Adata(immune_sub)
# Save the AnnData object to a local file
immune_adata_path <- file.path("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta", "immune_adata.h5ad")
adata.Save(immune_adata_path)

# Get cell counts per cell type per tissue
cell_counts <- table(placenta_seurat$celltype19, placenta_seurat$tissue_type)
cell_counts <- cell_counts[sort(rownames(cell_counts)),]
# cell_counts <- cell_counts[rownames(cell_counts)!=c("B cell-1", "B cell-2", "CD4 T-1", "CD4 T-2"),]
cell_counts <- as.data.frame(cell_counts)
cell_counts <- rownames_to_column(cell_counts, "cluster")

ggplot(data=cell_counts, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill") + 
  labs(title = "Cell Proportion Per Cluster Per Tissue Type", x= "Cluster", y= "Cell Proportion") + 
  guides(fill = guide_legend(title = "Tissue Type",reverse=FALSE)) +
  scale_fill_discrete(labels = c("DISC", "DEC", "CAM")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ddx3y <- Percent_Expressing(placenta_seurat, features = "DDX3Y")
ddx3y <- t(ddx3y)
ddx3y <- as.data.frame(ddx3y)
ddx3y$cluster <- levels(placenta_seurat)
ggplot(data=ddx3y, aes(x=cluster, y=DDX3Y)) +
  geom_bar(stat = "identity", fill = '#619CFF') + 
  labs(title = "Percent of Cells Expressing DDX3Y Per Cluster", x= "Cluster", y= "Percent Expressing") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Macro DDX3Y plot for paper
FeaturePlot_scCustom(macro_sub, features = "DDX3Y")

# Main UMAP for the paper used this code to generate it:
DimPlot_scCustom(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, colors_use = DiscretePalette_scCustomize(num_colors = 9,
                                                                                                                                                   palette = "varibow", shuffle_pal = FALSE), color_seed = 5)
