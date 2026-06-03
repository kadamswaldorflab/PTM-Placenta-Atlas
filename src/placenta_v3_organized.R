# =============================================================================
# Placenta Atlas v3 — Full Analysis Pipeline
# Join parts: cat placenta_v3_part*.R > placenta_v3_organized.R
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------

setwd("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta/placenta_v3")

library(Seurat)
library(DropletUtils)
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(scCustomize)
library(harmony)
library(future)
library(ggplot2)
library(purrr)
library(SeuratExtend)
library(SeuratDisk)
library(dittoSeq)
library(FNN)
library(ggrepel)
library(shadowtext)
library(gatepoints)
library(shiny)
library(plotly)
library(sccomp)
library(dplyr)
library(patchwork)

set.seed(42)
gradient <- c("blue", "white", "red")
options(future.globals.maxSize = 1e20)

# -----------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

save_umap <- function(plot, filename, width = 15) {
  plot <- plot +
    coord_fixed(ratio = 0.75) +
    theme(aspect.ratio = 0.75)
  ggsave(filename, plot = plot, width = width, height = width * 0.75, device = "pdf")
}

save_subset_for_scvelo <- function(seurat_obj, name) {
  metadata <- seurat_obj@meta.data
  metadata$barcode <- rownames(metadata)
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  umap_coords$barcode <- rownames(umap_coords)
  combined <- merge(metadata, umap_coords, by = "barcode")
  write.csv(combined, paste0(name, "_metadata.csv"), row.names = FALSE)
  cat("Saved metadata for", name, "with", nrow(combined), "cells\n")
}

qc_summary <- function(obj, stage) {
  md <- obj@meta.data
  data.frame(
    Stage         = stage,
    N_cells       = nrow(md),
    Mean_nFeature = round(mean(md$nFeature_RNA, na.rm = TRUE), 1),
    SD_nFeature   = round(sd(md$nFeature_RNA,   na.rm = TRUE), 1),
    Mean_nCount   = round(mean(md$nCount_RNA,   na.rm = TRUE), 1),
    SD_nCount     = round(sd(md$nCount_RNA,     na.rm = TRUE), 1),
    Mean_pct_mt   = round(mean(md$percent.mt,   na.rm = TRUE), 2),
    SD_pct_mt     = round(sd(md$percent.mt,     na.rm = TRUE), 2)
  )
}

# Interactive Plotly lasso gating — saves selected barcodes to gated_cells.rds
launch_lasso_app <- function(seurat_obj) {
  umap_coords    <- as.data.frame(Embeddings(seurat_obj, "umap"))
  cluster_labels <- as.character(Idents(seurat_obj))

  ui <- fluidPage(
    plotlyOutput("umap", height = "600px"),
    verbatimTextOutput("selected"),
    actionButton("save", "Save Selected Cells")
  )

  server <- function(input, output, session) {
    output$umap <- renderPlotly({
      plot_ly(umap_coords, x = ~umap_1, y = ~umap_2,
              type = "scattergl", mode = "markers",
              marker = list(size = 3),
              text = cluster_labels,
              hovertemplate = "Cluster: %{text}<extra></extra>",
              source = "umap") %>%
        layout(dragmode = "lasso")
    })
    selected_cells <- reactive({ event_data("plotly_selected", source = "umap") })
    output$selected <- renderPrint({
      req(selected_cells())
      paste("Selected:", nrow(selected_cells()), "cells")
    })
    observeEvent(input$save, {
      req(selected_cells())
      cells <- rownames(umap_coords)[selected_cells()$pointNumber + 1]
      saveRDS(cells, "gated_cells.rds")
      showNotification("Cells saved to gated_cells.rds")
    })
  }
  shinyApp(ui, server)
}

# Contour UMAP colored by tissue with cluster labels
plot_contour_umap <- function(seurat_obj, cell_type_col, padding = 1.5, top_padding = NULL,
                               title = "Subset UMAP with cell type labels",
                               filename = NULL, width = 15) {
  tissue_colors <- c("CV" = "#E41A1C", "Decidua" = "#4DAF4A", "MemInoc" = "#377EB8")
  tissue_labels <- c("CV" = "Disc",    "Decidua" = "Dec",     "MemInoc" = "Mem")
  y_top <- if (is.null(top_padding)) padding else top_padding

  umap_df <- data.frame(
    UMAP_1  = Embeddings(seurat_obj, reduction = "umap")[, 1],
    UMAP_2  = Embeddings(seurat_obj, reduction = "umap")[, 2],
    cluster = seurat_obj[[cell_type_col]][, 1],
    tissue  = seurat_obj$tissue
  )
  label_coords <- umap_df %>%
    group_by(cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  x_range <- range(umap_df$UMAP_1)
  y_range <- range(umap_df$UMAP_2)

  p <- DimPlot(seurat_obj, group.by = "tissue", reduction = "umap",
               cols = tissue_colors, pt.size = 0.5) +
    scale_color_manual(values = tissue_colors, labels = tissue_labels) +
    geom_density_2d(data = umap_df,
                    aes(x = UMAP_1, y = UMAP_2, group = cluster),
                    color = "black", linewidth = 0.5,
                    contour_var = "ndensity", breaks = 0.2, adjust = 2) +
    geom_text_repel(data = label_coords,
                    aes(x = UMAP_1, y = UMAP_2, label = cluster),
                    colour = "black", bg.color = "white", bg.r = 0.15,
                    size = 4, fontface = "bold", max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.3,
                    segment.color = "grey50", segment.size = 0.3,
                    min.segment.length = 0.2) +
    scale_x_continuous(limits = c(x_range[1] - padding, x_range[2] + padding)) +
    scale_y_continuous(limits = c(y_range[1] - padding, y_range[2] + y_top)) +
    coord_fixed(ratio = 0.75, clip = "off") +
    ggtitle(title) +
    theme(legend.position = "right", plot.margin = margin(10, 20, 10, 10)) +
    xlab("UMAP 1") + ylab("UMAP 2")

  if (!is.null(filename)) ggsave(filename, plot = p, width = width, height = width * 0.75)
  invisible(p)
}

# Polychrome UMAP with repelled shadow labels
plot_polychrome_umap <- function(seurat_obj, n_colors = 33, shuffle = TRUE,
                                  filename = NULL, width = 15, height = 10) {
  p <- DimPlot_scCustom(
    seurat_obj, reduction = "umap", label = FALSE, pt.size = 0.5, repel = TRUE,
    colors_use = DiscretePalette_scCustomize(
      num_colors = n_colors, palette = "polychrome", shuffle_pal = shuffle),
    color_seed = 42
  )
  label_data <- p$data %>%
    group_by(ident) %>%
    summarise(
      x = median(.data[[grep("_1$", colnames(p$data), value = TRUE)[1]]]),
      y = median(.data[[grep("_2$", colnames(p$data), value = TRUE)[1]]])
    )
  p2 <- p + geom_text_repel(
    data = label_data, aes(x = x, y = y, label = ident),
    colour = "black", bg.color = "white", bg.r = 0.15,
    size = 4, fontface = "bold", max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey50", segment.size = 0.3, min.segment.length = 0.2
  ) + xlab("UMAP 1") + ylab("UMAP 2")

  if (!is.null(filename)) ggsave(filename, plot = p2, width = width, height = height, device = "pdf")
  invisible(p2)
}

# Transfer subset labels back to the full object
transfer_labels_to_full <- function(full_obj, subset_obj, col) {
  full_obj[[col]] <- as.character(full_obj[[col]])
  full_obj[[col]][colnames(subset_obj)] <- as.character(Idents(subset_obj))
  full_obj[[col]] <- as.factor(full_obj[[col]])
  full_obj
}

# Reassign scattered cluster cells to nearest non-scattered neighbor by UMAP KNN
merge_scattered_cluster <- function(seurat_obj, scattered_cluster) {
  umap_coords  <- as.data.frame(Embeddings(seurat_obj, "umap"))
  cluster_ids  <- as.character(Idents(seurat_obj))
  names(cluster_ids) <- rownames(umap_coords)

  scattered_cells  <- rownames(umap_coords)[cluster_ids == scattered_cluster]
  reference_cells  <- rownames(umap_coords)[cluster_ids != scattered_cluster]
  nn <- get.knnx(data = umap_coords[reference_cells, ], query = umap_coords[scattered_cells, ], k = 1)

  new_labels <- cluster_ids
  for (i in seq_along(scattered_cells))
    new_labels[scattered_cells[i]] <- cluster_ids[reference_cells[nn$nn.index[i]]]

  Idents(seurat_obj) <- new_labels
  seurat_obj
}

# -----------------------------------------------------------------------------
# 3. LOAD DATA
# -----------------------------------------------------------------------------

placenta_seurat <- readRDS("placenta_seurat_FINAL_v4.RDS")
mac_subset      <- readRDS("mac_subset_v5.RDS")
t_subset        <- readRDS("t_subset_v3.RDS")
tb_subset       <- readRDS("tb_subset_FINAL.RDS")
dsc_subset      <- readRDS("dsc_subset.RDS")

# -----------------------------------------------------------------------------
# 4. PER-SAMPLE PROCESSING (emptyDrops → SoupX → QC)
# -----------------------------------------------------------------------------

data_dir       <- "/mmfs1/gscratch/kawaldorflab/jcorn427/placenta/placenta_v2/cell_ranger_outs"
sample_folders <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
data_dirs      <- setNames(
  lapply(sample_folders, function(x) file.path(x, "outs")),
  basename(sample_folders)
)

placenta_seuratect_list <- list()
for (sample_name in names(data_dirs)) {
  cat("Processing sample:", sample_name, "\n")

  sample_path   <- data_dirs[[sample_name]]
  raw_path      <- file.path(sample_path, "raw_feature_bc_matrix")
  filtered_path <- file.path(sample_path, "filtered_feature_bc_matrix")

  if (!dir.exists(raw_path) || !dir.exists(filtered_path)) {
    cat("Warning: Required matrix paths not found for", sample_name, "- skipping\n\n")
    next
  }

  tryCatch({
    cat("  Loading counts...\n")
    raw_counts      <- Seurat::Read10X(raw_path)
    filtered_counts <- Seurat::Read10X(filtered_path)

    # emptyDrops filtering
    cat("  Running emptyDrops...\n")
    set.seed(42)
    e.out   <- DropletUtils::emptyDrops(raw_counts)
    is_cell <- which(e.out$FDR < 0.01)
    cat("  emptyDrops called", length(is_cell), "cells (vs", ncol(filtered_counts), "from CellRanger)\n")

    if (length(is_cell) < 100) {
      cat("  Too few cells from emptyDrops, falling back to CellRanger filtered matrix\n")
    } else {
      filtered_counts <- raw_counts[, is_cell]
    }

    if (ncol(filtered_counts) == 0 || nrow(filtered_counts) == 0) {
      cat("Warning: No cells or features for", sample_name, "- skipping\n\n")
      next
    }

    cat("  Sample has", ncol(filtered_counts), "cells and", nrow(filtered_counts), "features\n")

    use_soupx     <- ncol(filtered_counts) >= 100
    soupx_success <- FALSE
    corrected_counts <- NULL

    if (use_soupx) {
      cat("  Creating preliminary Seurat object...\n")
      placenta_temp <- CreateSeuratObject(counts = filtered_counts, min.cells = 1, min.features = 1)
      placenta_temp <- NormalizeData(placenta_temp, verbose = FALSE)
      placenta_temp <- FindVariableFeatures(placenta_temp, verbose = FALSE)
      placenta_temp <- ScaleData(placenta_temp, verbose = FALSE)
      placenta_temp <- RunPCA(placenta_temp, verbose = FALSE)
      placenta_temp <- FindNeighbors(placenta_temp, dims = 1:10, verbose = FALSE)
      placenta_temp <- FindClusters(placenta_temp, resolution = 0.5, verbose = FALSE)
      prelim_clusters <- placenta_temp$seurat_clusters

      cat("  Attempting SoupX correction...\n")
      tryCatch({
        sc <- SoupChannel(raw_counts, filtered_counts)
        sc <- setClusters(sc, prelim_clusters)
        sc <- autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.9)
        corrected_counts <- adjustCounts(sc, roundToInt = TRUE)
        soupx_success <- TRUE
        cat("  SoupX correction successful!\n")
      }, error = function(e) {
        cat("  SoupX failed:", conditionMessage(e), "\n  Proceeding without SoupX...\n")
      })
    } else {
      cat("  Too few cells for SoupX (", ncol(filtered_counts), "). Skipping.\n")
    }

    fb <- CreateSeuratObject(
      counts  = if (soupx_success) corrected_counts else filtered_counts,
      project = sample_name
    )
    fb$sample          <- sample_name
    fb$soupx_corrected <- soupx_success

    fb[["percent.mt"]] <- PercentageFeatureSet(fb, pattern = "^MT-")
    cells_before <- ncol(fb)
    fb <- subset(fb, subset = nFeature_RNA > 200 & percent.mt < 5)
    cells_after  <- ncol(fb)

    if (cells_after == 0) {
      cat("Warning: No cells passed QC for", sample_name, "- skipping\n\n")
      next
    }

    placenta_seuratect_list[[sample_name]] <- fb
    cat("Completed:", sample_name, "with", cells_after, "cells (",
        cells_before - cells_after, "filtered out)",
        ifelse(soupx_success, "(SoupX corrected)", "(no SoupX)"), "\n\n")

  }, error = function(e) {
    cat("ERROR processing", sample_name, ":\n  ", conditionMessage(e), "\n\n")
  })
}

cat("\n========================================\n")
cat("Successfully processed", length(placenta_seuratect_list), "out of", length(data_dirs), "samples\n")

if (length(placenta_seuratect_list) > 0) {
  soupx_status <- sapply(placenta_seuratect_list, function(x) x$soupx_corrected[1])
  cat("SoupX corrected:", sum(soupx_status), "samples\n")
  cat("Not corrected:", sum(!soupx_status), "samples\n\n")
  cat("Cell counts per sample:\n")
  for (sname in names(placenta_seuratect_list))
    cat("  ", sname, ":", ncol(placenta_seuratect_list[[sname]]), "cells\n")
}

# -----------------------------------------------------------------------------
# 5. SAMPLE ID STANDARDIZATION & MERGE
# -----------------------------------------------------------------------------

sample_ids <- sapply(names(placenta_seuratect_list), function(x) {
  animal_id <- regmatches(x, regexpr("[JLRZAM]\\d+", x))
  tissue <- NA
  if      (grepl("_CV_|_cv_", x))                  tissue <- "CV"
  else if (grepl("de[cd]idua|De[cd]idua", x))      tissue <- "Decidua"
  else if (grepl("[Mm]em[Ii]noc", x))               tissue <- "MemInoc"
  if (length(animal_id) > 0 && !is.na(tissue)) paste0(animal_id, "_", tissue) else x
})

for (i in seq_along(placenta_seuratect_list))
  placenta_seuratect_list[[i]]@meta.data$sample <- sample_ids[i]

# Remove problematic sample
placenta_seuratect_list <- placenta_seuratect_list[!grepl("R06147", names(placenta_seuratect_list))]
sample_ids <- names(placenta_seuratect_list)

placenta_merged <- merge(
  x = placenta_seuratect_list[[1]],
  y = placenta_seuratect_list[2:length(placenta_seuratect_list)],
  add.cell.ids = sample_ids,
  project = "all_samples"
)

# -----------------------------------------------------------------------------
# 6. DOUBLET DETECTION & QC FILTERING
# -----------------------------------------------------------------------------

placenta_merged <- JoinLayers(placenta_merged)
sce <- as.SingleCellExperiment(placenta_merged)
sce <- scDblFinder(sce, samples = "orig.ident")
placenta_merged$scDblFinder_class <- sce$scDblFinder.class
placenta_merged$scDblFinder_score <- sce$scDblFinder.score
placenta_merged <- PercentageFeatureSet(placenta_merged, pattern = "^MT", col.name = "percent.mt")

VlnPlot(placenta_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample")
ggsave("qc_plots_pre_filt.pdf", width = 30, height = 10)

plot1 <- FeatureScatter(placenta_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(placenta_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

placenta_merged_filt <- subset(placenta_merged,
                               subset = nFeature_RNA > 200 &
                                        nFeature_RNA < 7500 &
                                        percent.mt  < 2.5  &
                                        scDblFinder_class == "singlet")

VlnPlot(placenta_merged_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample")

# Rename barcodes to sample_barcode format
new_barcodes <- paste0(placenta_merged_filt$sample, "_",
                       sub(".*_([ACGT]+-[0-9]+)$", "\\1", colnames(placenta_merged_filt)))
colnames(placenta_merged_filt) <- new_barcodes
rownames(placenta_merged_filt@meta.data) <- new_barcodes

placenta_merged_filt$tissue <- sapply(strsplit(placenta_merged_filt$sample, "_"), function(x) x[2])

head(colnames(placenta_merged_filt))
head(placenta_merged_filt@meta.data[, c("sample", "tissue")])

# -----------------------------------------------------------------------------
# 7. NORMALIZATION, PCA & HARMONY INTEGRATION
# -----------------------------------------------------------------------------

placenta_merged_filt <- SCTransform(placenta_merged_filt, vars.to.regress = "percent.mt")
placenta_merged_filt <- RunPCA(placenta_merged_filt, assay = "SCT", npcs = 50)
ElbowPlot(placenta_merged_filt, ndims = 50, reduction = "pca")

DimPlot(placenta_merged_filt, reduction = "pca", group.by = "sample")
ggsave("pca_pre_harmony_sample.pdf", width = 15, height = 10)

placenta_merged_filt <- RunUMAP(placenta_merged_filt, reduction = "pca", dims = 1:30,
                                reduction.name = "umap_preharmony",
                                reduction.key  = "UMAPpreharmony_")
p <- DimPlot(placenta_merged_filt, reduction = "umap_preharmony", group.by = "sample", shuffle = TRUE)
save_umap(p, "umap_preharmony.pdf")

saveRDS(placenta_merged_filt, file = "placenta_merged_pre_harmony.RDS")

placenta_merged_filt <- RunHarmony(
  object         = placenta_merged_filt,
  group.by.vars  = "sample",
  reduction.use  = "pca",
  dims.use       = 1:30,
  theta          = 2,
  max.iter.harmony = 20,
  verbose        = TRUE
)

placenta_merged_filt <- RunUMAP(placenta_merged_filt, reduction = "harmony", dims = 1:30,
                                n.neighbors = 30, min.dist = 0.3)
placenta_merged_filt <- FindNeighbors(placenta_merged_filt, reduction = "harmony", dims = 1:30)
placenta_merged_filt <- FindClusters(placenta_merged_filt, resolution = 0.5)

p1 <- DimPlot(placenta_merged_filt, reduction = "umap", group.by = "sample") +
  ggtitle("By Sample - After Harmony")
save_umap(p1, "post_harmony_umap.pdf")

DimPlot(placenta_merged_filt, reduction = "harmony", group.by = "sample")
ggsave("pca_post_harmony_sample.pdf", width = 15, height = 10)

# QC metrics summary table
qc_table <- bind_rows(
  qc_summary(placenta_merged,      "Pre-filtering"),
  qc_summary(placenta_merged_filt, "Post-filtering & Harmonization"),
  qc_summary(placenta_seurat,      "Final (annotated)")
)
print(qc_table)
write.csv(qc_table, "qc_metrics_table.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. MACROPHAGE SUBSET ANNOTATION
# -----------------------------------------------------------------------------

macro_features <- c("CEACAM8", "SELL", "CSF3R", "CD177", "OLFM4", "CLEC7A", "ADAMDEC1", "IL12B", "THBS1", "MAMU-DRB1",
                    "CD80", "IL10", "IL1B", "IL6", "TNF", "CD209", "SPP1", "FOLR2", "FCGR1A", "HMOX1",
                    "LYVE1", "CCL2", "DDX3Y", "RPS4Y1", "MRC1", "CD163", "CD68", "CD14", "CCL17",
                    "CCR2", "S100A8", "S100A9", "LYZ", "MS4A7", "ITGAM")

clusters_of_interest <- c("HB", "FMAC", "MAC1", "MAC2", "Mo-MAC1", "Mo-MAC2", "NEUT", "MAC3")
mac_subset <- subset(placenta_seurat, idents = clusters_of_interest)
mac_subset <- SCTransform(mac_subset)
mac_subset <- RunPCA(mac_subset, npcs = 30)
ElbowPlot(mac_subset, ndims = 30)
mac_subset <- RunHarmony(mac_subset, group.by.vars = "sample", reduction = "pca",
                          assay.use = "SCT", reduction.save = "harmony")
mac_subset <- FindNeighbors(mac_subset, reduction = "harmony", dims = 1:30)
mac_subset <- FindClusters(mac_subset, resolution = 0.15)
mac_subset <- RunUMAP(mac_subset, reduction = "harmony", dims = 1:30)

p1 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE) +
  ggtitle("Macrophage Subset — Resolution 0.15")
ggsave("mac_umap_1.pdf", width = 15, height = 10, units = "in", dpi = 300)

macro_B <- DotPlot_scCustom(mac_subset, features = macro_features,
                             remove_axis_titles = FALSE, flip_axes = TRUE,
                             x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
macro_B_filtered <- macro_B +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
macro_B_filtered

p1 | macro_B_filtered
ggsave("macro_plots.pdf", width = 25, height = 10)

p2 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v4") +
  ggtitle("Main UMAP calls")
p1 | p2
ggsave("macro_plots_2.pdf", width = 25, height = 10)

DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "tissue")

# --- Round 1: lasso gate to label HB cells ---
launch_lasso_app(mac_subset)
gated_cells <- readRDS("gated_cells.rds")
mac_subset$updated_clusters <- as.character(Idents(mac_subset))
mac_subset$updated_clusters[gated_cells] <- "HB"
mac_subset$updated_clusters <- as.factor(mac_subset$updated_clusters)
Idents(mac_subset) <- "updated_clusters"
table(mac_subset$updated_clusters)

mac_subset <- RenameIdents(mac_subset,
                            "0" = "M2", "1" = "MDM", "3" = "M1", "4" = "NEUT",
                            "6" = "M2a MAC", "7" = "FMAC", "9" = "M2")
mac_subset$cell_type_v5 <- Idents(mac_subset)
DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v5")

# --- Round 2: relabel cluster 2 via lasso ---
launch_lasso_app(mac_subset)
gated_cells <- readRDS("gated_cells.rds")
mac_subset$cell_type_v5 <- as.character(mac_subset$cell_type_v5)
mac_subset$cell_type_v5[gated_cells] <- "2"
mac_subset$cell_type_v5 <- as.factor(mac_subset$cell_type_v5)
Idents(mac_subset) <- "cell_type_v5"
table(mac_subset$cell_type_v5)

# --- Merge scattered cluster 5 to nearest neighbor ---
mac_subset <- merge_scattered_cluster(mac_subset, "5")
mac_subset$cell_type_v5 <- Idents(mac_subset)
saveRDS(mac_subset, "mac_subset_merged.rds")

# --- Sub-cluster refinement: cluster 8 ---
mac_subset$cell_type_v4 <- Idents(mac_subset)
mac_subset <- FindSubCluster(mac_subset, "8", graph.name = "SCT_nn", resolution = 0.25)

p1 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "sub.cluster") +
  ggtitle("Macrophage Subset — 8 Sub")
ggsave("mac_umap_1.pdf", width = 15, height = 10, units = "in", dpi = 300, plot = p1)

macro_B <- DotPlot_scCustom(mac_subset, features = macro_features, group.by = "sub.cluster",
                             remove_axis_titles = FALSE, flip_axes = TRUE,
                             x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
macro_B_filtered <- macro_B +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
macro_B_filtered

# --- Sub-cluster cluster 3 ---
mac_subset$cell_type_v4 <- mac_subset$sub.cluster
Idents(mac_subset) <- "cell_type_v4"
mac_subset <- FindSubCluster(mac_subset, "3", graph.name = "SCT_nn", resolution = 0.2)

p1 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "sub.cluster") +
  ggtitle("Macrophage Subset — 8, 3 Sub")
ggsave("mac_umap_1.pdf", width = 15, height = 10, units = "in", dpi = 300, plot = p1)

macro_B <- DotPlot_scCustom(mac_subset, features = macro_features, group.by = "sub.cluster",
                             remove_axis_titles = FALSE, flip_axes = TRUE,
                             x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
macro_B_filtered <- macro_B +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
macro_B_filtered
p1 | macro_B_filtered
ggsave("macro_plots.pdf", width = 25, height = 10)

Idents(mac_subset) <- "sub.cluster"
mac_subset <- RenameIdents(mac_subset, "1" = "PAMM", "10" = "M2a", "11" = "Mono")
mac_subset <- RenameIdents(mac_subset, "3_3" = "3_0", "3_4" = "3_0", "8_4" = "8_3")
mac_subset <- RenameIdents(mac_subset, "Mono" = "11")
mac_subset <- RenameIdents(mac_subset, "3_1" = "Mono")
mac_subset <- RenameIdents(mac_subset, "3_2" = "NEUT", "3_0" = "NEUT")
mac_subset <- RenameIdents(mac_subset, "2" = "0")
Cluster_Highlight_Plot(mac_subset, "5", reduction = "umap")

# --- Sub-cluster cluster 5 ---
mac_subset$cell_type_v5 <- Idents(mac_subset)
Idents(mac_subset) <- "cell_type_v5"
mac_subset <- FindSubCluster(mac_subset, "5", graph.name = "SCT_nn", resolution = 0.05)

p1 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE, group.by = "sub.cluster") +
  ggtitle("Macrophage Subset — 8, 3, 5 Sub")
ggsave("mac_umap_1.pdf", width = 15, height = 10, units = "in", dpi = 300, plot = p1)

Idents(mac_subset) <- "sub.cluster"
mac_subset <- RenameIdents(mac_subset, "9" = "ADAMDEC1+MAC")
Cluster_Highlight_Plot(mac_subset, "3_5", reduction = "umap")

p1 <- DimPlot_scCustom(mac_subset, reduction = "umap", label = TRUE)
ggsave("mac_umap_1.pdf", width = 15, height = 10, units = "in", dpi = 300, plot = p1)

macro_B <- DotPlot_scCustom(mac_subset, features = macro_features,
                             remove_axis_titles = FALSE, flip_axes = TRUE,
                             x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
macro_B_filtered <- macro_B +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
macro_B_filtered

mac_subset <- RenameIdents(mac_subset, "5_0" = "SPP1+MAC")
mac_subset <- RenameIdents(mac_subset, "4" = "M1", "5_1" = "0")
mac_subset <- RenameIdents(mac_subset, "8_0" = "FMAC", "8_1" = "HB")
mac_subset <- RenameIdents(mac_subset, "8_3" = "LYVE+MAC", "8_2" = "LYVE+MAC")
mac_subset <- RenameIdents(mac_subset, "3_5" = "PAMM", "5_2" = "0", "5_3" = "0")

# --- Marker analysis ---
mac_subset <- PrepSCTFindMarkers(mac_subset)
all_markers <- FindAllMarkers(object = mac_subset) %>% Add_Pct_Diff()

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = mac_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = mac_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 13)

top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_immune.csv")

pdf("mac_mark.pdf", width = 20, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = mac_subset, features = top_markers,
                               flip = TRUE, x_lab_rotate = 90, k = 13)
print(dot_plot)
dev.off()

saveRDS(mac_subset, file = "mac_subset_v2.RDS")

# --- Transfer mac labels to full object ---
placenta_seurat <- transfer_labels_to_full(placenta_seurat, mac_subset, "cell_type_v4")
Idents(placenta_seurat) <- "cell_type_v4"
table(placenta_seurat$cell_type_v4)
length(colnames(mac_subset)) == sum(colnames(placenta_seurat) %in% colnames(mac_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v4")
ggsave("umap_v9.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Marker analysis on full object
placenta_seurat <- PrepSCTFindMarkers(placenta_seurat)
all_markers <- FindAllMarkers(object = placenta_seurat) %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 27)
top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_immune.csv")

pdf("placenta_clust.pdf", width = 30, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                               flip = TRUE, x_lab_rotate = 90, k = 27)
print(dot_plot)
dev.off()

clusters_of_interest <- c("0", "11", "6", "7")
top_markers <- all_markers %>%
  filter(cluster %in% clusters_of_interest) %>%
  Extract_Top_Markers(marker_dataframe = ., num_genes = 7, named_vector = FALSE,
                      make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 5)

pdf("placenta_clust_unk.pdf", width = 15, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                               flip = TRUE, x_lab_rotate = 90, k = 5)
print(dot_plot)
dev.off()

saveRDS(mac_subset, file = "mac_subset_v3.RDS")
saveRDS(placenta_seurat, file = "placenta_mapped_v4.RDS")

# -----------------------------------------------------------------------------
# 11. T / NK / B CELL SUBSET ANNOTATION
# -----------------------------------------------------------------------------

t_features <- c("CD3D", "CD4", "ENSMMUG00000003532", "CD8A", "FCGR3",
                 "NCAM1", "TRAC", "ICOS", "CD69", "ITK",
                 "THEMIS", "GRAP2", "CD28", "CD6", "BCL11B",
                 "ITGAE", "IL7R", "PDCD1", "IFNG", "FOXP3",
                 "IL2RA", "ITGA1", "NKG2A", "CD94", "KLRB1",
                 "NKG7", "NKG2D", "PRF1", "GZMA", "GZMB",
                 "GZMK", "GZMM", "TIGIT", "CCL5", "ITGA4",
                 "CD38", "CD7", "KLRD1", "ZBTB38", "ZBTB16")

Idents(t_subset) <- "cell_type_v4"

Fig_7D <- DotPlot_scCustom(t_subset, features = t_features, group.by = "cell_type_v4",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_7D_filtered <- Fig_7D +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_7D_filtered
ggsave("t_dot.pdf", width = 10, height = 10)

# Lasso gate unknown T cells
launch_lasso_app(t_subset)
gated_cells <- readRDS("gated_cells.rds")
t_subset$cell_type_v4 <- as.character(t_subset$cell_type_v4)
t_subset$cell_type_v4[gated_cells] <- "unknown"
t_subset$cell_type_v4 <- as.factor(t_subset$cell_type_v4)
Idents(t_subset) <- "cell_type_v4"
table(t_subset$cell_type_v4)

t_subset <- PrepSCTFindMarkers(t_subset)
all_markers <- FindAllMarkers(object = t_subset) %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = t_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = t_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 27)
top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_immune.csv")

top_markers_unk <- Extract_Top_Markers(
  marker_dataframe = all_markers %>% filter(cluster == "unknown"),
  num_genes = 50, named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC"
)

DimPlot_scCustom(t_subset, reduction = "umap", group.by = "tissue")

# B cell dot plot (full object)
b_features <- c("CD38", "MS4A1", "CD19", "CD79A", "PAX5",
                 "EBF1", "BLK", "BANK1", "BACH2", "AFF3",
                 "CLEC17A", "IGHM", "ENSMMUG00000015202", "JCHAIN", "ENSMMUG00000059850",
                 "ENSMMUG00000044861", "ENSMMUG00000040771", "ENSMMUG00000002764", "TNFRSF17", "MZB1")
Fig_7D <- DotPlot_scCustom(placenta_seurat, features = b_features, group.by = "cell_type_v4",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_7D_filtered <- Fig_7D +
  scale_y_discrete(limits = c("unknown", "Plasma")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_7D_filtered
ggsave("b_dot.pdf", width = 10, height = 10)

# Transfer T subset labels to full object
placenta_seurat <- transfer_labels_to_full(placenta_seurat, t_subset, "cell_type_v4")
Idents(placenta_seurat) <- "cell_type_v4"
table(placenta_seurat$cell_type_v4)
length(colnames(t_subset)) == sum(colnames(placenta_seurat) %in% colnames(t_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v4")
ggsave("umap_v11.pdf", width = 15, height = 10, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# 12. TROPHOBLAST SUBSET ANNOTATION
# -----------------------------------------------------------------------------

DimPlot_scCustom(tb_subset, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "tissue")

# Lasso gate unknown trophoblasts
launch_lasso_app(tb_subset)
gated_cells <- readRDS("gated_cells.rds")
tb_subset$cell_type_v4 <- as.character(tb_subset$cell_type_v4)
tb_subset$cell_type_v4[gated_cells] <- "unknown_tb"
tb_subset$cell_type_v4 <- as.factor(tb_subset$cell_type_v4)
Idents(tb_subset) <- "cell_type_v4"
table(tb_subset$cell_type_v4)

# Transfer TB labels to full object
placenta_seurat <- transfer_labels_to_full(placenta_seurat, tb_subset, "cell_type_v4")
Idents(placenta_seurat) <- "cell_type_v4"
table(placenta_seurat$cell_type_v4)
length(colnames(tb_subset)) == sum(colnames(placenta_seurat) %in% colnames(tb_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v4")
ggsave("umap_v12.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Trophoblast marker dot plot
tb_features <- c("TFAP2C", "GATA3", "ITGA6", "EPCAM", "CDH1",
                  "SPINT1", "TP63", "PEG3", "KRT5", "KRT7",
                  "TFAP2A", "SLC2A1", "IGF2", "KRT19", "SLC1A5",
                  "GCM1", "YAP1", "PAGE4", "CDH5", "CYP19A1",
                  "CGA", "CGB2", "CGB1–8", "CSH1", "GH2",
                  "PSGs", "SLC38A", "CK6C", "KLF4", "KRT17",
                  "PEG10", "CD56", "CD31", "CD146", "KRT8",
                  "KRT18", "EpCAM", "PCAM1", "CSH2", "CYP11A1",
                  "CRH", "PAPPA", "PAX8", "HNF1B", "PAEP", "EHF", "GRHL2", "DPP4", "AGR2",
                  "LEPR", "CLDN3", "CLDN4", "SCGB1A1", "RXFP1", "MUC1", "MUC16", "GLYCAM1", "GPX3")

Fig_4D <- DotPlot_scCustom(placenta_seurat, features = tb_features, group.by = "cell_type_v4",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_4D_filtered <- Fig_4D +
  scale_y_discrete(limits = c("schCTB", "CTB1", "CTB2", "CTB3", "STB", "unknown_tb")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_4D_filtered
ggsave("4d_dot_tb.pdf", width = 8, height = 10)

markers <- FindMarkers(object = placenta_seurat, ident.1 = "unknown_tb",
                       min.pct = 0.25, logfc.threshold = 0.25)
top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  pull(gene)
markers <- rownames_to_column(markers, var = "gene")
write_csv(markers, file = "unknown_tb_markers.csv")

FeaturePlot_scCustom(placenta_seurat, features = "LEPR", reduction = "umap")

# Rename and save v5
placenta_seurat <- RenameIdents(placenta_seurat,
                                "CD4 T" = "Treg", "NK" = "pNK", "7" = "M0", "0" = "M2",
                                "6" = "0", "11" = "MKI67+MAC",
                                "unknown" = "B", "unknown_tb" = "GE", "Plasma" = "PB")
placenta_seurat$cell_type_v5 <- Idents(placenta_seurat)

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v13.pdf", width = 15, height = 10, units = "in", dpi = 300)
saveRDS(placenta_seurat, file = "placenta_mapped_v5.RDS")

mac_subset <- RenameIdents(mac_subset, "7" = "M0", "0" = "M2", "6" = "0", "11" = "MKI67+MAC")
mac_subset$cell_type_v5 <- Idents(mac_subset)
saveRDS(mac_subset, file = "mac_subset_v3.RDS")

t_subset <- RenameIdents(t_subset, "unknown" = "B", "CD4 T" = "Treg", "NK" = "pNK")
t_subset$cell_type_v5 <- Idents(t_subset)
saveRDS(t_subset, file = "t_subset_v3.RDS")

tb_subset <- RenameIdents(tb_subset, "unknown_tb" = "GE")
tb_subset$cell_type_v5 <- Idents(tb_subset)
saveRDS(tb_subset, file = "tb_subset_v3.RDS")

# Re-cluster trophoblast subset (CTB / STB)
clusters_of_interest <- c("schCTB", "CTB2", "CTB1", "CTB3", "STB")
tb_subset <- subset(placenta_seurat, idents = clusters_of_interest)
tb_subset <- SCTransform(tb_subset)
tb_subset <- RunPCA(tb_subset, npcs = 30)
ElbowPlot(tb_subset, ndims = 30)
tb_subset <- RunHarmony(tb_subset, group.by.vars = "sample", reduction = "pca",
                         assay.use = "SCT", reduction.save = "harmony")
tb_subset <- FindNeighbors(tb_subset, reduction = "harmony", dims = 1:30)
tb_subset <- FindClusters(tb_subset, resolution = 0.15)
tb_subset <- RunUMAP(tb_subset, reduction = "harmony", dims = 1:30)

p1 <- DimPlot_scCustom(tb_subset, reduction = "umap", label = TRUE, group.by = "tissue")
ggsave("tb_umap_2.pdf", plot = p1, width = 15, height = 10, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# 13. IMMUNE SUBSET RE-CLUSTERING (NK / T / B)
# -----------------------------------------------------------------------------

imm_features <- c("CD3D", "CD4", "ENSMMUG0000003532", "CD8A",
                   "FCGR3", "NCAM1", "ITGAM", "ITGA1", "ENTPD1",
                   "ENSMMUG00000013289", "KLRB1", "PRF1", "ENSMMUG00000063583", "GZMA",
                   "GZMB", "NKG7", "CTSW", "GZMM", "XCL1",
                   "ITGAE", "BICDL1", "STAT4", "NKG2D", "CD57",
                   "CSFR2", "EOMES", "CD39", "CD73", "ADORA3",
                   "CSF1", "CSF1R", "CCR1", "XCR1",
                   "TRAC", "ICOS", "CD69", "ITK",
                   "THEMIS", "GRAP2", "CD28", "CD6", "BCL11B",
                   "IL7R", "PDCD1", "IFNG", "FOXP3",
                   "IL2RA", "NKG2A", "CD94",
                   "GZMK", "TIGIT", "CCL5", "ITGA4",
                   "CD38", "CD7", "KLRD1", "ZBTB38", "ZBTB16",
                   "MZB1", "TNFRSF17", "ENSMMUG00000002764", "ENSMMUG00000040771", "ENSMMUG00000044861",
                   "ENSMMUG00000059850", "JCHAIN", "ENSMMUG00000015202", "IGHM", "CLEC17A",
                   "AFF3", "BACH2", "BANK1", "BLK", "EBF1",
                   "PAX5", "CD79A", "CD19", "MS4A1", "HAVCR2", "LAG3", "CTLA4", "CD244", "TNFRSF9",
                   "TOX", "TOX2", "BATF", "IRF4", "NR4A2", "NR4A3", "IKZF2", "IKZF3",
                   "PD-1", "TCF7", "TCF1", "MS4A2", "CPA3", "FCER1A")

clusters_of_interest <- c("B", "pNK", "uNK1", "uNK2", "uNK3", "Treg", "TRM", "CD8 TRM", "PB")
imm_subset <- subset(placenta_seurat, idents = clusters_of_interest)
imm_subset <- SCTransform(imm_subset)
imm_subset <- RunPCA(imm_subset, npcs = 30)
ElbowPlot(imm_subset, ndims = 30)
imm_subset <- RunHarmony(imm_subset, group.by.vars = "sample", reduction = "pca",
                          assay.use = "SCT", reduction.save = "harmony")
imm_subset <- FindNeighbors(imm_subset, reduction = "harmony", dims = 1:30)
imm_subset <- FindClusters(imm_subset, resolution = 0.15)
imm_subset <- RunUMAP(imm_subset, reduction = "harmony", dims = 1:30)

p1 <- DimPlot_scCustom(imm_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v5")
ggsave("imm_umap_1.pdf", plot = p1, width = 15, height = 10, units = "in", dpi = 300)

# Misc stromal/endothelial dot plot (Fig 9)
misc_features_v1 <- c("PDGFRB", "MCAM", "MYH11", "ACTA2", "THY1",
                       "COL6A2", "DCN", "LUM", "VIM", "FN1",
                       "PECAM1", "VWF", "CDH5", "TEK", "KDR",
                       "PAX8", "HNF1B", "PAEP", "SCGB1A1", "EHF",
                       "CDH1", "EPCAM", "CLDN3", "CLDN4", "GRHL2", "KRT19", "AGR3", "RXFP1")
Fig_9 <- DotPlot_scCustom(placenta_seurat, features = misc_features_v1, group.by = "cell_type_v5",
                           remove_axis_titles = FALSE, flip_axes = TRUE,
                           x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_9_filtered <- Fig_9 +
  scale_y_discrete(limits = c("PC", "FIB", "EC", "GE")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_9_filtered
ggsave("misc_dot.pdf", plot = Fig_9_filtered, width = 10, height = 10)

# Lasso gate unknown immune cells
launch_lasso_app(imm_subset)
gated_cells <- readRDS("gated_cells.rds")
imm_subset$cell_type_v5 <- as.character(imm_subset$cell_type_v5)
imm_subset$cell_type_v5[gated_cells] <- "unknown"
imm_subset$cell_type_v5 <- as.factor(imm_subset$cell_type_v5)
Idents(imm_subset) <- "cell_type_v5"
table(imm_subset$cell_type_v5)

imm_subset <- PrepSCTFindMarkers(imm_subset)
all_markers <- FindAllMarkers(object = imm_subset) %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = imm_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = imm_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 11)
top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_immune.csv")

pdf("imm_mark.pdf", width = 20, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = imm_subset, features = top_markers,
                               flip = TRUE, x_lab_rotate = 90, k = 11)
print(dot_plot)
dev.off()

unk_markers <- all_markers %>% filter(cluster == "unknown") %>% arrange(desc(avg_log2FC))
write.csv(unk_markers, file = "unknown_imm.csv")

FeaturePlot(imm_subset, features = c("MS4A2", "CPA3", "CD6", "ZAP70"), ncol = 2, reduction = "umap")
FeatureScatter(imm_subset, feature1 = "MS4A2", feature2 = "CD6")

# Transfer immune labels to full object (round 1)
placenta_seurat <- transfer_labels_to_full(placenta_seurat, imm_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(imm_subset)) == sum(colnames(placenta_seurat) %in% colnames(imm_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v13.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Sub-cluster unknown cells
imm_subset <- FindSubCluster(imm_subset, cluster = "unknown", graph.name = "SCT_snn", resolution = 0.05)
DimPlot_scCustom(imm_subset, group.by = "sub.cluster", label = TRUE, reduction = "umap")
ggsave("immune_umap.pdf", width = 15, height = 10)

Idents(imm_subset) <- imm_subset$sub.cluster
placenta_seurat <- transfer_labels_to_full(placenta_seurat, imm_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(imm_subset)) == sum(colnames(placenta_seurat) %in% colnames(imm_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v14.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Marker analysis for unknown sub-clusters
for (unk_id in c("unknown_0", "unknown_1", "unknown_2")) {
  markers_unk <- FindMarkers(object = placenta_seurat, ident.1 = unk_id,
                             min.pct = 0.25, logfc.threshold = 0.25)
  markers_unk <- rownames_to_column(markers_unk, var = "gene")
  write_csv(markers_unk, file = paste0(sub("_", "", unk_id), "_markers.csv"))
}

# NK/T/B dot plot with unknown subclusters
Fig_6D <- DotPlot_scCustom(placenta_seurat, features = imm_features, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_6D_filtered <- Fig_6D +
  scale_y_discrete(limits = c("pNK", "uNK1", "uNK2", "uNK3", "TRM", "CD8 TRM", "Treg",
                               "B", "PB", "unknown_0", "unknown_1", "unknown_2")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_6D_filtered
ggsave("dot_nk_t_b.pdf", width = 10, height = 16)

# Further sub-cluster unknown_1
imm_subset <- FindSubCluster(imm_subset, cluster = "unknown_1", graph.name = "SCT_snn", resolution = 0.5)
DimPlot_scCustom(imm_subset, group.by = "sub.cluster", label = TRUE, reduction = "umap")
ggsave("immune_umap.pdf", width = 15, height = 10)
table(imm_subset$sub.cluster)

Fig_6D <- DotPlot_scCustom(imm_subset, features = imm_features, group.by = "sub.cluster",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_6D_filtered <- Fig_6D +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_6D_filtered
ggsave("dot_nk_t_b.pdf", width = 10, height = 16)

Idents(imm_subset) <- imm_subset$sub.cluster
placenta_seurat <- transfer_labels_to_full(placenta_seurat, imm_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(imm_subset)) == sum(colnames(placenta_seurat) %in% colnames(imm_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v14.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Expanded dot plot with unknown_1 subclusters
imm_features_ext <- c(imm_features, "MKI67", "FLT3", "CXCR1", "CD34")
Fig_6D <- DotPlot_scCustom(placenta_seurat, features = imm_features_ext, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_6D_filtered <- Fig_6D +
  scale_y_discrete(limits = c("pNK", "uNK1", "uNK2", "uNK3", "TRM", "CD8 TRM", "Treg",
                               "B", "PB", "unknown_0", "unknown_1_0", "unknown_1_1",
                               "unknown_1_2", "unknown_1_3", "unknown_2")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_6D_filtered
ggsave("dot_nk_t_b.pdf", width = 10, height = 16)

# Collapse unknowns to ImmP
imm_subset <- RenameIdents(imm_subset,
                            "unknown_2"   = "ImmP", "unknown_0"   = "ImmP",
                            "unknown_1_3" = "ImmP", "unknown_1_0" = "ImmP",
                            "unknown_1_2" = "ImmP", "unknown_1_1" = "ImmP",
                            "CD8 TRM"     = "CD8+TRM")
imm_subset$cell_type_v5 <- Idents(imm_subset)

placenta_seurat <- transfer_labels_to_full(placenta_seurat, imm_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(imm_subset)) == sum(colnames(placenta_seurat) %in% colnames(imm_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v15.pdf", width = 15, height = 10, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# 14. FINAL RENAMING & LABEL PROPAGATION
# -----------------------------------------------------------------------------

placenta_seurat <- RenameIdents(placenta_seurat, "LYVE+MAC" = "LYVE1+MAC")
placenta_seurat$cell_type_v5 <- Idents(placenta_seurat)
saveRDS(placenta_seurat, file = "placenta_seurat_FINAL.RDS")

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "tissue")
ggsave("umap_v15.pdf", width = 15, height = 10, units = "in", dpi = 300)

mac_subset <- RenameIdents(mac_subset, "LYVE+MAC" = "LYVE1+MAC", "0" = "M2")
mac_subset$cell_type_v5 <- Idents(mac_subset)
saveRDS(mac_subset, file = "mac_subset_v4.RDS")

mac_subset <- RenameIdents(mac_subset, "PAMM" = "MDM")
mac_subset$cell_type_v5 <- Idents(mac_subset)
mac_subset <- RenameIdents(mac_subset, "M2" = "M-INT")
mac_subset$cell_type_v5 <- Idents(mac_subset)

placenta_seurat <- transfer_labels_to_full(placenta_seurat, mac_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(mac_subset)) == sum(colnames(placenta_seurat) %in% colnames(mac_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v16.pdf", width = 15, height = 10, units = "in", dpi = 300)

saveRDS(placenta_seurat, "placenta_seurat_FINAL_v2.RDS")
saveRDS(mac_subset,      "mac_subset_v5.RDS")
saveRDS(tb_subset,       "tb_subset_v4.RDS")

# Add fetal sex metadata
Idents(placenta_seurat) <- "id"
placenta_seurat <- RenameIdents(placenta_seurat,
                                "CTRL01" = "Female", "CTRL02" = "Male",  "CTRL03" = "Female",
                                "CTRL04" = "Male",   "CTRL05" = "Male",  "CTRL07" = "Male",
                                "CTRL08" = "Female", "CTRL09" = "Female","CTRL10" = "Female")
placenta_seurat$fet_sex <- Idents(placenta_seurat)

# -----------------------------------------------------------------------------
# 15. TROPHOBLAST FINAL RE-CLUSTERING & CLEANUP
# -----------------------------------------------------------------------------

# Lasso gate additional unknown trophoblasts
Idents(tb_subset) <- "cell_type_v5"
launch_lasso_app(tb_subset)
gated_cells <- readRDS("gated_cells.rds")
tb_subset$cell_type_v5 <- as.character(tb_subset$cell_type_v5)
tb_subset$cell_type_v5[gated_cells] <- "unknown4"
tb_subset$cell_type_v5 <- as.factor(tb_subset$cell_type_v5)
Idents(tb_subset) <- "cell_type_v5"
table(tb_subset$cell_type_v5)

DimPlot_scCustom(tb_subset, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_tb.pdf", width = 15, height = 10, units = "in", dpi = 300)

tb_subset <- PrepSCTFindMarkers(tb_subset)
all_markers <- FindAllMarkers(object = tb_subset) %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = tb_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = tb_subset, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 11)
top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_tb.csv")

pdf("tb_mark.pdf", width = 20, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = tb_subset, features = top_markers,
                               flip = TRUE, x_lab_rotate = 90, k = 10)
print(dot_plot)
dev.off()

placenta_seurat <- transfer_labels_to_full(placenta_seurat, tb_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)
length(colnames(tb_subset)) == sum(colnames(placenta_seurat) %in% colnames(tb_subset))

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v17.pdf", width = 15, height = 10, units = "in", dpi = 300)

# Full object marker analysis
placenta_seurat <- PrepSCTFindMarkers(placenta_seurat)
all_markers <- FindAllMarkers(object = placenta_seurat) %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7,
                                   named_vector = FALSE, make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = placenta_seurat, features = top_markers,
                            flip = TRUE, x_lab_rotate = 90, k = 27)
top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20,
                                      data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers_tb.csv")

# Top markers for unresolved TB populations
clusters_of_interest <- c("CTB1", "CTB2", "CTB3", "schCTB", "unknown1", "unknown2", "unknown3", "unknown4")
gene_cluster_map <- all_markers %>%
  filter(cluster %in% clusters_of_interest) %>%
  group_by(gene) %>%
  summarise(other_clusters = paste(cluster, collapse = ", "), .groups = "drop")
top_markers <- all_markers %>%
  filter(cluster %in% clusters_of_interest) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  distinct(gene, .keep_all = TRUE) %>%
  ungroup() %>%
  left_join(gene_cluster_map, by = "gene") %>%
  mutate(other_clusters = mapply(function(oc, cl) {
    paste(setdiff(strsplit(oc, ", ")[[1]], cl), collapse = ", ")
  }, other_clusters, cluster))
write.csv(top_markers, "top_markers.csv", row.names = FALSE)

FeaturePlot_scCustom(tb_subset, features = "FOXO1", reduction = "umap")

# DSC dot plot
dsc_features <- c("CCND1", "PDGFRA", "PDGFRB", "HOXA11", "POSTN",
                   "TIMP3", "IGFBP1", "RGCC", "CXCL12", "APOD",
                   "PGR", "HAND2", "WT1", "FOXO1", "COL4A5",
                   "DKK1", "PRL", "WNT5A", "COL1A1", "COL3A1",
                   "VIM", "THY1", "CD90", "BMP2", "WNT4", "LEFTY2")
Fig_8D <- DotPlot_scCustom(placenta_seurat, features = dsc_features, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_8D_filtered <- Fig_8D +
  scale_y_discrete(limits = c("ESC", "DSC1", "DSC2", "unknown1", "unknown2", "unknown3", "unknown4")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_8D_filtered
ggsave("8d_dot_dsc.pdf", width = 8, height = 10)

# Relabel unknown TB by KNN majority vote
nn_graph    <- tb_subset@graphs$SCT_snn
cell_labels <- tb_subset@meta.data$cell_type_v5
names(cell_labels) <- colnames(tb_subset)

schCTB_cells <- which(cell_labels == "schCTB")
new_labels   <- cell_labels
for (cell_idx in schCTB_cells) {
  neighbors       <- which(nn_graph[cell_idx, ] > 0)
  neighbor_labels <- cell_labels[neighbors]
  neighbor_labels <- neighbor_labels[neighbor_labels != "schCTB"]
  if (length(neighbor_labels) > 0)
    new_labels[cell_idx] <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
}
tb_subset@meta.data$cell_type_relabeled <- new_labels
tb_subset$cell_type_relabeled <- droplevels(tb_subset$cell_type_relabeled)
table(tb_subset$cell_type_relabeled)

# Inspect mast cell markers in unknown4
VlnPlot(tb_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        idents = c("unknown4", "CTB1", "CTB2", "CTB3", "STB"),
        group.by = "cell_type_v5", pt.size = 0.1)

all_markers %>%
  filter(cluster == "unknown4") %>%
  arrange(desc(avg_log2FC)) %>%
  head(30) %>%
  print(n = 30)

FeaturePlot_scCustom(tb_subset,
                     features = c("MS4A2", "CPA3", "FCER1A", "KIT", "ENSMMUG00000015103"),
                     reduction = "umap", num_columns = 3)
ggsave("mastmarkers.pdf", width = 15, height = 10)

# Remove mast cell contamination (unknown4)
unknown4_cells  <- colnames(tb_subset)[tb_subset$cell_type_v5 == "unknown4"]
tb_subset       <- subset(tb_subset, cell_type_v5 != "unknown4")
tb_subset$cell_type_v5 <- droplevels(tb_subset$cell_type_v5)

cells_to_keep   <- colnames(placenta_seurat)[!colnames(placenta_seurat) %in% unknown4_cells]
placenta_seurat <- subset(placenta_seurat, cells = cells_to_keep)
cat("Removed", length(unknown4_cells), "cells\n")
cat("placenta_seurat now has", ncol(placenta_seurat), "cells\n")

table(tb_subset$cell_type_v5)
table(tb_subset$cell_type_relabeled)
tb_subset$cell_type_relabeled <- droplevels(tb_subset$cell_type_relabeled)

Idents(tb_subset) <- "cell_type_relabeled"
tb_subset <- RenameIdents(tb_subset, "unknown1" = "AEC", "unknown2" = "FIB2", "unknown3" = "schCTB")
tb_subset$cell_type_v5 <- Idents(tb_subset)
table(tb_subset$cell_type_v5)

placenta_seurat <- transfer_labels_to_full(placenta_seurat, tb_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v18.pdf", width = 15, height = 10, units = "in", dpi = 300)

saveRDS(tb_subset,       "tb_subset_v5.RDS")
saveRDS(placenta_seurat, "placenta_seurat_FINAL_v3.RDS")

# Swap CTB1 <-> CTB3 labels in both objects
for (obj_name in c("placenta_seurat", "tb_subset")) {
  obj <- get(obj_name)
  obj$cell_type_v5 <- as.character(obj$cell_type_v5)
  obj$cell_type_v5[obj$cell_type_v5 == "CTB1"]     <- "CTB_TEMP"
  obj$cell_type_v5[obj$cell_type_v5 == "CTB3"]     <- "CTB1"
  obj$cell_type_v5[obj$cell_type_v5 == "CTB_TEMP"] <- "CTB3"
  obj$cell_type_v5 <- as.factor(obj$cell_type_v5)
  Idents(obj) <- "cell_type_v5"
  assign(obj_name, obj)
}
table(placenta_seurat$cell_type_v5)
table(tb_subset$cell_type_v5)

# Rename FIB -> FIB1 in full object
placenta_seurat$cell_type_v5 <- as.character(placenta_seurat$cell_type_v5)
placenta_seurat$cell_type_v5[placenta_seurat$cell_type_v5 == "FIB"] <- "FIB1"
placenta_seurat$cell_type_v5 <- as.factor(placenta_seurat$cell_type_v5)
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)

# Final TB re-clustering
tb_clusters <- c("CTB1", "CTB2", "CTB3", "STB", "schCTB")
tb_subset <- subset(placenta_seurat, idents = tb_clusters)
table(tb_subset$cell_type_v5)

tb_subset <- SCTransform(tb_subset)
tb_subset <- RunPCA(tb_subset, npcs = 30)
tb_subset <- RunHarmony(tb_subset, group.by.vars = "sample", reduction = "pca",
                         assay.use = "SCT", reduction.save = "harmony")
tb_subset <- FindNeighbors(tb_subset, reduction = "harmony", dims = 1:30)
tb_subset <- FindClusters(tb_subset, resolution = 0.15)
tb_subset <- RunUMAP(tb_subset, reduction = "harmony", dims = 1:30)

DimPlot_scCustom(tb_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v5")
ggsave("tb_umap_final.pdf", width = 15, height = 10)

saveRDS(tb_subset, "tb_subset_FINAL.RDS")

# Lasso gate final CTB1 corrections
launch_lasso_app(tb_subset)
gated_cells <- readRDS("gated_cells.rds")
tb_subset$cell_type_v5 <- as.character(tb_subset$cell_type_v5)
tb_subset$cell_type_v5[gated_cells] <- "CTB1"
tb_subset$cell_type_v5 <- as.factor(tb_subset$cell_type_v5)
Idents(tb_subset) <- "cell_type_v5"
table(tb_subset$cell_type_v5)

DimPlot_scCustom(tb_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v5")

placenta_seurat <- transfer_labels_to_full(placenta_seurat, tb_subset, "cell_type_v5")
Idents(placenta_seurat) <- "cell_type_v5"
table(placenta_seurat$cell_type_v5)

DimPlot_scCustom(placenta_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")
ggsave("umap_v19.pdf", width = 15, height = 10, units = "in", dpi = 300)

saveRDS(tb_subset,       "tb_subset_FINAL.RDS")
saveRDS(placenta_seurat, "placenta_seurat_FINAL_v4.RDS")

# Trajectory subset (CTB → STB)
traj_clusters <- c("CTB1", "CTB2", "CTB3", "STB")
traj_subset <- subset(placenta_seurat, idents = traj_clusters)
traj_subset <- SCTransform(traj_subset)
traj_subset <- RunPCA(traj_subset, npcs = 30)
ElbowPlot(traj_subset, ndims = 30)
traj_subset <- RunHarmony(traj_subset, group.by.vars = "sample", reduction = "pca",
                           assay.use = "SCT", reduction.save = "harmony")
traj_subset <- FindNeighbors(traj_subset, reduction = "harmony", dims = 1:30)
traj_subset <- FindClusters(traj_subset, resolution = 0.15)
traj_subset <- RunUMAP(traj_subset, reduction = "harmony", dims = 1:30)

DimPlot_scCustom(traj_subset, reduction = "umap", label = TRUE, group.by = "cell_type_v5")
ggsave("traj_umap_v1.pdf", width = 15, height = 10)

# Export scVelo metadata
save_subset_for_scvelo(mac_subset, "mac")
save_subset_for_scvelo(tb_subset,  "tb")

metadata    <- placenta_seurat@meta.data
metadata$barcode <- rownames(metadata)
umap_coords <- as.data.frame(placenta_seurat@reductions$umap@cell.embeddings)
umap_coords$barcode <- rownames(umap_coords)
combined    <- merge(metadata, umap_coords, by = "barcode")
write.csv(combined, "full_placenta_metadata.csv", row.names = FALSE)
cat("Saved metadata for full dataset with", nrow(combined), "cells\n")

# -----------------------------------------------------------------------------
# 16. TREATMENT & SAMPLE METADATA
# -----------------------------------------------------------------------------

treatment_df <- data.frame(
  animal_id = c("Z15006", "A10007", "M10190", "Z12353", "L10152",
                 "Z15060", "Z15198", "J08186", "R10195"),
  treatment = c("Media",  "Media",  "Media",  "Media",  "Media",
                 "Saline", "Saline", "Saline", "Saline")
)

meta <- placenta_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(animal_id = sub("_.*", "", sample)) %>%
  left_join(treatment_df, by = "animal_id")
placenta_seurat$animal_id <- meta$animal_id
placenta_seurat$treatment  <- meta$treatment
table(placenta_seurat$treatment, placenta_seurat$animal_id)

# -----------------------------------------------------------------------------
# 17. PAPER FIGURES
# -----------------------------------------------------------------------------

Idents(placenta_seurat) <- "cell_type_v5"

tissue_colors <- c("CV" = "#E41A1C", "Decidua" = "#4DAF4A", "MemInoc" = "#377EB8")
tissue_labels <- c("CV" = "Disc",    "Decidua" = "Dec",     "MemInoc" = "Mem")

# --- Full atlas: tissue UMAP without contours ---
umap_df_full <- data.frame(
  UMAP_1  = Embeddings(placenta_seurat, reduction = "umap")[, 1],
  UMAP_2  = Embeddings(placenta_seurat, reduction = "umap")[, 2],
  cluster = placenta_seurat$cell_type_v5,
  tissue  = placenta_seurat$tissue
)
label_coords_full <- umap_df_full %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
x_range <- range(umap_df_full$UMAP_1)
y_range <- range(umap_df_full$UMAP_2)
padding <- 1.5

p_full <- DimPlot(placenta_seurat, group.by = "tissue", reduction = "umap",
                  cols = tissue_colors, pt.size = 0.5) +
  scale_color_manual(values = tissue_colors, labels = tissue_labels) +
  geom_text_repel(data = label_coords_full,
                  aes(x = UMAP_1, y = UMAP_2, label = cluster),
                  colour = "black", bg.color = "white", bg.r = 0.15,
                  size = 4, fontface = "bold", max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey50", segment.size = 0.3, min.segment.length = 0.2) +
  scale_x_continuous(limits = c(x_range[1] - padding, x_range[2] + padding)) +
  scale_y_continuous(limits = c(y_range[1] - padding, y_range[2] + padding)) +
  coord_fixed(ratio = 0.75, clip = "off") +
  ggtitle("Tissue UMAP with cell type labels") +
  theme(legend.position = "right", plot.margin = margin(10, 20, 10, 10)) +
  xlab("UMAP 1") + ylab("UMAP 2")
w <- 15
ggsave("no_contour_plot.pdf", plot = p_full, width = w, height = w * 0.75)

# --- Full atlas polychrome UMAP ---
plot_polychrome_umap(placenta_seurat, n_colors = 38, shuffle = FALSE, filename = "polychrome_true.pdf")

# --- Subset contour UMAPs (tissue-colored, density contours per cluster) ---
plot_contour_umap(mac_subset,  "cell_type_v5",              filename = "contour_umap_mac.pdf")
plot_contour_umap(tb_subset,   "cell_type_v5", top_padding = 3.5, filename = "contour_umap_tb.pdf")
plot_contour_umap(imm_subset,  "cell_type_v5",              filename = "contour_umap_imm.pdf")
plot_contour_umap(dsc_subset,  "cell_type_v4",              filename = "contour_umap_dsc.pdf")

# --- Subset polychrome UMAPs ---
plot_polychrome_umap(mac_subset, filename = "umap_plot_mac.pdf")
plot_polychrome_umap(tb_subset,  filename = "umap_plot_tb_weird.pdf")
plot_polychrome_umap(imm_subset, filename = "umap_plot_imm.pdf")

# --- RPS4Y1 expression bar (all cells) ---
plot_df <- data.frame(
  cell_type  = placenta_seurat$cell_type_v5,
  expression = GetAssayData(placenta_seurat, layer = "data")["RPS4Y1", ]
) %>%
  group_by(cell_type) %>%
  summarise(pct_expressing = mean(expression > 0) * 100) %>%
  mutate(cell_type = factor(cell_type, levels = sort(unique(as.character(cell_type)))))

ggplot(plot_df, aes(x = cell_type, y = pct_expressing)) +
  geom_bar(stat = "identity", fill = "#2166AC") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Cell Type", y = "% Cells Expressing RPS4Y1", title = "RPS4Y1 Expression")
ggsave("RPS4Y1_percent_bar.pdf", height = 10, width = 15)

# --- Tissue proportion bar ---
plot_df <- placenta_seurat@meta.data %>%
  group_by(cell_type_v5, tissue) %>%
  summarise(n = n()) %>%
  group_by(cell_type_v5) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(cell_type_v5 = factor(cell_type_v5, levels = sort(unique(as.character(cell_type_v5)))))

ggplot(plot_df, aes(x = cell_type_v5, y = proportion, fill = tissue)) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(labels = c("Decidua" = "Dec", "CV" = "Disc", "MemInoc" = "Mem")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Proportion", fill = "Tissue Type")
ggsave("tissue_percent_bar.pdf", height = 10, width = 15)

# --- DDX3Y bar (male cells only) ---
male_seurat <- subset(placenta_seurat, subset = fet_sex == "Male")
plot_df <- data.frame(
  cell_type  = male_seurat$cell_type_v5,
  expression = GetAssayData(male_seurat, layer = "data")["DDX3Y", ]
) %>%
  group_by(cell_type) %>%
  summarise(pct_expressing = mean(expression > 0) * 100) %>%
  mutate(cell_type = factor(cell_type, levels = sort(unique(as.character(cell_type)))))

ggplot(plot_df, aes(x = cell_type, y = pct_expressing)) +
  geom_bar(stat = "identity", fill = "#2166AC") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Cell Type", y = "% Cells Expressing DDX3Y", title = "DDX3Y Expression (Male)")
ggsave("DDX3Y_percent_bar_male.pdf", height = 10, width = 15)

# --- NK / T / B dot plot — Fig 6D ---
imm_features_fig <- c("CD3D", "CD4", "ENSMMUG0000003532", "CD8A",
                       "FCGR3", "NCAM1", "ITGAM", "ITGA1", "ENTPD1",
                       "ENSMMUG00000013289", "KLRB1", "PRF1", "ENSMMUG00000063583", "GZMA",
                       "GZMB", "NKG7", "CTSW", "GZMM", "XCL1",
                       "ITGAE", "BICDL1", "STAT4", "NKG2D", "CD57",
                       "CSFR2", "EOMES", "CD39", "CD73", "ADORA3",
                       "CSF1", "CSF1R", "CCR1", "XCR1",
                       "TRAC", "ICOS", "CD69", "ITK",
                       "THEMIS", "GRAP2", "CD28", "CD6", "BCL11B",
                       "IL7R", "PDCD1", "IFNG", "FOXP3",
                       "IL2RA", "NKG2A", "CD94",
                       "GZMK", "TIGIT", "CCL5", "ITGA4",
                       "CD38", "CD7", "KLRD1", "ZBTB38", "ZBTB16",
                       "MZB1", "TNFRSF17", "ENSMMUG00000002764", "ENSMMUG00000040771", "ENSMMUG00000044861",
                       "ENSMMUG00000059850", "JCHAIN", "ENSMMUG00000015202", "IGHM", "CLEC17A",
                       "AFF3", "BACH2", "BANK1", "BLK", "EBF1",
                       "PAX5", "CD79A", "CD19", "MS4A1", "HAVCR2", "LAG3", "CTLA4", "CD244", "TNFRSF9",
                       "TOX", "TOX2", "BATF", "IRF4", "NR4A2", "NR4A3", "IKZF2", "IKZF3",
                       "PD-1", "TCF7", "TCF1", "MS4A2", "CPA3", "FCER1A", "MKI67", "FLT3", "CXCR1", "CD34")

Fig_6D <- DotPlot_scCustom(placenta_seurat, features = imm_features_fig, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_6D_filtered <- Fig_6D +
  scale_y_discrete(limits = c("pNK", "uNK1", "uNK2", "uNK3", "TRM", "CD8 TRM", "Treg", "B", "PB", "ImmP")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_6D_filtered
ggsave("dot_nk_t_b.pdf", width = 10, height = 16)

# --- Macrophage dot plot — Fig 5D ---
mac_features_fig <- c("S100A8", "S100A9", "VCAN", "LST1", "SELL",
                       "CCR2", "THBS1", "IL1B", "CD62L",
                       "MAMU-DRB1", "CD86", "CD9", "IL10", "CD209",
                       "FOLR2", "FCGR1A", "HMOX1", "CCL2",
                       "CLEC7A", "MRC1", "CD163", "CD68", "CD14",
                       "CSF1R", "C1QB", "CCL17", "ADAMDEC1", "MKI67", "LYVE1",
                       "SPP1", "DDX3Y", "RPS4Y1", "ARG1", "CCL1")
Fig_5D <- DotPlot_scCustom(placenta_seurat, features = mac_features_fig, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_5D_filtered <- Fig_5D +
  scale_y_discrete(limits = c("NEUT", "Mono", "MDM", "M0", "M1", "M-INT", "M2a",
                               "ADAMDEC1+MAC", "MKI67+MAC", "LYVE1+MAC", "SPP1+MAC", "FMAC", "HB")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_5D_filtered
ggsave("5d_dot_mac.pdf", width = 8, height = 10)

# --- Trophoblast dot plot — Fig 4D ---
tb_features_fig <- c("TFAP2C", "GATA3", "ITGA6", "EPCAM", "CDH1",
                      "SPINT1", "TP63", "PEG3", "KRT5", "KRT7",
                      "TFAP2A", "SLC2A1", "IGF2", "KRT19", "SLC1A5",
                      "GCM1", "YAP1", "PAGE4", "CDH5", "CYP19A1",
                      "CGA", "CGB2", "CGB1–8", "CSH1", "GH2",
                      "PSGs", "SLC38A", "CK6C", "KLF4", "KRT17",
                      "PEG10", "CD56", "CD31", "CD146", "KRT8",
                      "KRT18", "EpCAM", "PCAM1", "CSH2", "CYP11A1",
                      "CRH", "PAPPA", "VTCN1", "LGALS16", "KRT18", "KRT6A")
Fig_4D <- DotPlot_scCustom(placenta_seurat, features = tb_features_fig, group.by = "cell_type_v5",
                            remove_axis_titles = FALSE, flip_axes = TRUE,
                            x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_4D_filtered <- Fig_4D +
  scale_y_discrete(limits = c("schCTB", "CTB1", "CTB2", "CTB3", "STB",
                               "unknown1", "unknown2", "unknown3", "unknown4")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_4D_filtered
ggsave("4d_dot_tb.pdf", width = 8, height = 10)

# --- Stromal/endothelial dot plot — Fig 9 ---
misc_features <- c("CDH1", "EPCAM", "CLDN4", "KRT19", "AGR2", "AGR3", "KRT24", "KRT16", "MSLN", "LAMB3",
                    "VTCN1", "LGALSL", "OVOL2", "COL6A2", "COL3A1", "COL1A1", "PAX8", "HNF1B", "PAEP",
                    "SCGB1A1", "EHF", "RXFP1", "ITGB3", "DPP4", "MAL", "MAL3",
                    "PDGFRB", "MCAM", "MYH11", "ACTA2", "THY1",
                    "DCN", "LUM", "VIM", "FN1", "DKK1", "IGFBP1", "POSTN",
                    "HBA", "HBE1", "ALAS2", "PECAM1", "VWF", "CDH5", "TEK", "KDR")
Fig_9 <- DotPlot_scCustom(placenta_seurat, features = misc_features, group.by = "cell_type_v5",
                           remove_axis_titles = FALSE, flip_axes = TRUE,
                           x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")
Fig_9_filtered <- Fig_9 +
  scale_y_discrete(limits = c("AEC", "GE", "PC", "FIB1", "FIB2", "EB", "EC")) +
  guides(size  = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
Fig_9_filtered
ggsave("9_dot_misc.pdf", width = 8, height = 10)

# -----------------------------------------------------------------------------
# 18. COMPOSITIONAL ANALYSIS (scComp)
# -----------------------------------------------------------------------------

# Treatment: Media vs Saline
sccomp_result <- placenta_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  sccomp_estimate(formula_composition = ~ treatment, .sample = sample,
                  .cell_group = cell_type_v5, bimodal_mean_variability_association = TRUE, cores = 4)
sccomp_result <- sccomp_result %>% sccomp_test(contrasts = c("treatmentSaline"))
print(sccomp_result)

plot_1d <- sccomp_result %>% plot_1D_intervals()
ggsave("sccomp_1D_intervals.pdf", plot = plot_1d, width = 10, height = 8)
plot_box <- sccomp_result %>% sccomp_boxplot(factor = "treatment")
ggsave("sccomp_boxplot.pdf", plot = plot_box, width = 12, height = 10)

# Tissue comparison: CV vs Decidua vs MemInoc
table(placenta_seurat$tissue, placenta_seurat$cell_type_v5)
sccomp_result_tissue <- placenta_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  sccomp_estimate(formula_composition = ~ tissue, .sample = sample,
                  .cell_group = cell_type_v5, bimodal_mean_variability_association = TRUE, cores = 4)
sccomp_result_tissue <- sccomp_result_tissue %>%
  sccomp_test(contrasts = c("tissueDecidua", "tissueMemInoc"))
print(sccomp_result_tissue)

plot_1d_tissue <- sccomp_result_tissue %>% plot_1D_intervals()
ggsave("sccomp_1D_intervals_tissue.pdf", plot = plot_1d_tissue, width = 10, height = 8)
plot_box_tissue <- sccomp_result_tissue %>% sccomp_boxplot(factor = "tissue")
ggsave("sccomp_boxplot_tissue.pdf", plot = plot_box_tissue, width = 12, height = 10)

# Fetal sex: Male vs Female
table(placenta_seurat$fet_sex, placenta_seurat$cell_type_v5)
sccomp_result_sex <- placenta_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  sccomp_estimate(formula_composition = ~ fet_sex, .sample = sample,
                  .cell_group = cell_type_v5, bimodal_mean_variability_association = TRUE, cores = 4)
sccomp_result_sex <- sccomp_result_sex %>% sccomp_test(contrasts = c("fet_sexMale"))
print(sccomp_result_sex)

plot_1d_sex <- sccomp_result_sex %>% plot_1D_intervals()
ggsave("sccomp_1D_intervals_sex.pdf", plot = plot_1d_sex, width = 10, height = 8)
plot_box_sex <- sccomp_result_sex %>% sccomp_boxplot(factor = "fet_sex")
ggsave("sccomp_boxplot_sex.pdf", plot = plot_box_sex, width = 12, height = 10)

# Export FDR-filtered results
cat("\n--- Treatment scComp: significant (c_FDR < 0.05) ---\n")
sccomp_result %>%
  select(cell_type_v5, parameter, c_lower, c_effect, c_upper, c_pH0, c_FDR) %>%
  filter(c_FDR < 0.05) %>% arrange(c_FDR) %>% print(n = Inf)
write.csv(
  sccomp_result %>% select(cell_type_v5, parameter, c_lower, c_effect, c_upper, c_pH0, c_FDR),
  "sccomp_treatment_FDR.csv", row.names = FALSE
)

cat("\n--- Fetal sex scComp: significant (c_FDR < 0.05) ---\n")
sccomp_result_sex %>%
  select(cell_type_v5, parameter, c_lower, c_effect, c_upper, c_pH0, c_FDR) %>%
  filter(c_FDR < 0.05) %>% arrange(c_FDR) %>% print(n = Inf)
write.csv(
  sccomp_result_sex %>% select(cell_type_v5, parameter, c_lower, c_effect, c_upper, c_pH0, c_FDR),
  "sccomp_sex_FDR.csv", row.names = FALSE
)
