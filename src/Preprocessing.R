library(dplyr)
library(igraph)
library(Seurat)
library(patchwork)
library(readr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(harmony)
library(SingleR)
library(pheatmap)
library(foreach)
library(doParallel)
library(distances)
library(scDEED)
library(celldex)

set.seed(42)

outFolder <- "./HumanNHP_mmulAnalysis/Final_Pipeline/"
dir.create(outFolder, showWarnings = FALSE, recursive = TRUE)

tissue_name <- "PLACENTA"

CtrlSal.data <- Read10X(data.dir = "/mmfs1/gscratch/kawaldorflab/rli/HumanNHP_mmulAnalysis/mmul_aggr/outs/count/filtered_feature_bc_matrix")
AggrCSV <- read.table("/mmfs1/gscratch/kawaldorflab/rli/HumanNHP_mmulAnalysis/mmul_aggr/outs/aggregation.csv", sep = ",", header = TRUE)
index <- as.numeric(sapply(strsplit(colnames(CtrlSal.data), "-"), function(X) X[2]))
AggrCSV <- AggrCSV[index, ]
rownames(AggrCSV) <- colnames(CtrlSal.data)

seu.obj <- CreateSeuratObject(counts = CtrlSal.data, project = paste0("CtrlSal", tissue_name), min.cells = 3, min.features = 200, meta.data = AggrCSV)
seu.obj$percent.mt <- PercentageFeatureSet(seu.obj, pattern = "^MT")

save(seu.obj, file = file.path(outFolder, "UnFilteredSeurat.Rdata"))

seu.filtered <- subset(
  seu.obj,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 800 & percent.mt < 2.5
)

qc_unf <- VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(filename = file.path(outFolder, paste0(tissue_name, "_qcMetricsGraph_Unfiltered.png")), plot = qc_unf, width = 20, height = 15, dpi = 300)

qc_filt <- VlnPlot(seu.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(filename = file.path(outFolder, paste0(tissue_name, "_qcMetricsGraph_Filtered.png")), plot = qc_filt, width = 20, height = 15, dpi = 300)

save(seu.filtered, file = file.path(outFolder, paste0(tissue_name, "FilteredSeurat_preHarmony.Rdata")))

seu.filtered <- SCTransform(seu.filtered, vars.to.regress = c("percent.mt"), verbose = FALSE)
seu.filtered <- RunPCA(seu.filtered, assay = "SCT", npcs = 30)

pca_elbow <- ElbowPlot(seu.filtered, ndims = 30, reduction = "pca")
ggsave(filename = file.path(outFolder, paste0(tissue_name, "_ElbowPlot_pca_30.png")), plot = pca_elbow, width = 12, height = 8, dpi = 300)

d <- 30
r <- 0.5

if (!"Subgroup" %in% colnames(seu.filtered@meta.data)) seu.filtered$Subgroup <- "UNK"
seu.filtered <- RunHarmony(seu.filtered, group.by.vars = "Subgroup")
seu.filtered <- FindNeighbors(seu.filtered, dims = 1:d, reduction = "harmony")
seu.filtered <- FindClusters(seu.filtered, resolution = r)
seu.filtered <- RunUMAP(seu.filtered, dims = 1:d, reduction = "harmony")

cond_plot <- DimPlot(seu.filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste0("UMAP by Clusters, dims=", d, ", res=", r))
ggsave(filename = file.path(outFolder, paste0(tissue_name, "_UMAP_byClusters_dims", d, "_res", r, ".png")), plot = cond_plot, width = 12, height = 10, dpi = 300)

treat_plot <- DimPlot(seu.filtered, reduction = "umap", split.by = "Subgroup", label = TRUE) +
  ggtitle(paste0("UMAP by Subgroup, dims=", d, ", res=", r))
ggsave(filename = file.path(outFolder, paste0(tissue_name, "_UMAP_bySubgroup_dims", d, "_res", r, ".png")), plot = treat_plot, width = 20, height = 10, dpi = 300)

markers <- FindAllMarkers(seu.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20)
write.csv(top20_markers, file.path(outFolder, paste0(tissue_name, "_top20_markers_per_cluster.csv")), row.names = FALSE)

save(seu.filtered, file = file.path(outFolder, paste0(tissue_name, "_FilteredSeurat_Harmony.Rdata")))

sl24 <- readr::read_rds("/mmfs1/gscratch/kawaldorflab/Globus/nardhy/third_tri_scRNA-Seq/sc.NormByLocationRep.Harmony.final.rds")

K <- 8
n_neighbors_parameters <- c(5, 20, 30, 40, 50)
min_dist_parameters <- c(0.1, 0.4)

Cell.Similarity <- function(pre_embedding_distances, pre_embedding_distances_permuted, reduced_dim_distances, reduced_dim_distances_permuted, similarity_percent = 0.5, makeCluster_n = 1) {
  numberselected <- floor((dim(pre_embedding_distances)[2]) * similarity_percent)
  cl <- parallel::makeCluster(makeCluster_n)
  doParallel::registerDoParallel(cl)
  rho_original <- foreach::foreach(i = 1:(dim(pre_embedding_distances)[2]), .combine = "c") %dopar% {
    require(distances)
    cor(
      (reduced_dim_distances[i, order(pre_embedding_distances[i, ])][2:(numberselected + 1)]),
      (sort(reduced_dim_distances[i, ])[2:(numberselected + 1)])
    )
  }
  rho_permuted <- foreach::foreach(i = 1:(dim(pre_embedding_distances)[2]), .combine = "c") %dopar% {
    require(distances)
    cor(
      (reduced_dim_distances_permuted[i, order(pre_embedding_distances_permuted[i, ])][2:(numberselected + 1)]),
      (sort(reduced_dim_distances_permuted[i, ])[2:(numberselected + 1)])
    )
  }
  parallel::stopCluster(cl)
  list(rho_original = as.numeric(rho_original), rho_permuted = as.numeric(rho_permuted))
}

input_data.permuted <- scDEED::Permuted(seu.filtered, K = K, slot = "scale.data", default_assay = "active.assay")
results.PCA <- scDEED::Distances.pre_embedding(seu.filtered, input_data.permuted, K = K, pre_embedding = "pca")

result_umap <- tibble()
ClassifiedCells <- NULL

for (n in n_neighbors_parameters) {
  for (m in min_dist_parameters) {
    results <- scDEED::Distances.UMAP(seu.filtered, input_data.permuted, K = K, pre_embedding = "pca", n = n, m = m, rerun = TRUE)
    similarity_score <- Cell.Similarity(
      results.PCA$pre_embedding_distances,
      results.PCA$pre_embedding_distances_permuted,
      results$reduced_dim_distances,
      results$reduced_dim_distances_permuted,
      similarity_percent = 0.5,
      makeCluster_n = 5
    )
    ClassifiedCells <- scDEED::Cell.Classify(similarity_score$rho_original, similarity_score$rho_permuted, dubious_cutoff = 0.05, trustworthy_cutoff = 0.95)
    dub <- ifelse(length(ClassifiedCells$dubious_cells) != 0, paste(ClassifiedCells$dubious_cells, collapse = ","), "none")
    int <- ifelse(length(ClassifiedCells$intermediate_cells) != 0, paste(ClassifiedCells$intermediate_cells, collapse = ","), "none")
    trust <- ifelse(length(ClassifiedCells$trustworthy_cells) != 0, paste(ClassifiedCells$trustworthy_cells, collapse = ","), "none")
    df <- tibble(n_neighbors = n, min.dist = m, number_dubious_cells = length(ClassifiedCells$dubious_cells), dubious_cells = dub, trustworthy_cells = trust, intermediate_cells = int)
    result_umap <- bind_rows(result_umap, df)
  }
}

readr::write_csv(result_umap, file.path(outFolder, "scDEED_param_sweep_results.csv"))
save.image(file = file.path(outFolder, "280225environment_woGraphs.RData"))

seu.filtered$classification <- "Unclassified"
if (!is.null(ClassifiedCells)) {
  if (length(ClassifiedCells$dubious_cells) > 0) seu.filtered$classification[ClassifiedCells$dubious_cells] <- "dubious"
  if (length(ClassifiedCells$trustworthy_cells) > 0) seu.filtered$classification[ClassifiedCells$trustworthy_cells] <- "trustworthy"
  if (length(ClassifiedCells$intermediate_cells) > 0) seu.filtered$classification[ClassifiedCells$intermediate_cells] <- "intermediate"
}

umap_plot <- DimPlot(
  seu.filtered, reduction = "umap", group.by = "classification",
  cols = c("dubious" = "red", "intermediate" = "grey", "trustworthy" = "blue", "Unclassified" = "black"), raster = FALSE
) +
  ggtitle("UMAP Plot of Cell Classifications") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_blank()) +
  labs(color = "Cell Classification")
ggsave(filename = file.path(outFolder, "Placenta_Dubious.png"), plot = umap_plot, width = 20, height = 15, dpi = 300)

idx <- which.min(result_umap$number_dubious_cells)
n_opt <- result_umap$n_neighbors[idx]
m_opt <- result_umap$min.dist[idx]

seu.filtered <- RunUMAP(seu.filtered, dims = 1:K, min.dist = m_opt, n.neighbors = n_opt, seed.use = 42)
opt_plot <- DimPlot(seu.filtered, reduction = "umap", group.by = "classification", cols = c("dubious" = "red", "intermediate" = "grey", "trustworthy" = "blue", "Unclassified" = "black"), raster = FALSE) +
  ggtitle("UMAP Plot of Cell Classifications (Optimized)")
ggsave(filename = file.path(outFolder, "Placenta_Dubious_optimized.png"), plot = opt_plot, width = 20, height = 15, dpi = 300)

cond_plot2 <- DimPlot(seu.filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP by Clusters")
ggsave(filename = file.path(outFolder, "Placenta_umap_unlabelled.png"), plot = cond_plot2, width = 20, height = 15, dpi = 300)

sl24_qc <- VlnPlot(sl24, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(filename = file.path(outFolder, "HumanPlacenta_qcMetricsGraph_Unfiltered.png"), plot = sl24_qc, width = 20, height = 15, dpi = 300)

sl24 <- subset(sl24, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 800 & percent.mt < 10)
sl24 <- SCTransform(sl24, vars.to.regress = c("percent.mt"), verbose = FALSE)
sl24 <- RunPCA(sl24, assay = "SCT", npcs = 100) %>% FindNeighbors(dims = 1:100) %>% FindClusters() %>% RunUMAP(dims = 1:100)

comb <- merge(x = seu.filtered, y = sl24)
comb <- SCTransform(comb, vars.to.regress = c("percent.mt"), verbose = FALSE)
comb <- ScaleData(comb)
comb <- RunPCA(comb, npcs = 100, verbose = FALSE) %>% FindNeighbors(dims = 1:100) %>% FindClusters() %>% RunUMAP(dims = 1:100)

reference <- celldex::HumanPrimaryCellAtlasData()
test_expr <- GetAssayData(seu.filtered, slot = "data")
pred <- SingleR(test = test_expr, ref = reference, labels = reference$label.main)

seu.filtered$SingleR.labels <- pred$labels[colnames(seu.filtered)]
anno_plot <- DimPlot(seu.filtered, reduction = "umap", group.by = "SingleR.labels", label = TRUE, raster = FALSE)
nonanno_plot <- DimPlot(seu.filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
ggsave(filename = file.path(outFolder, "Placenta_umap_labelled_SingleR.png"), plot = anno_plot, width = 20, hei





