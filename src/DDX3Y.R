library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

outFolder <- "./HumanNHP_mmulAnalysis/Final_Pipeline/final_figures"
dir.create(outFolder, showWarnings = FALSE, recursive = TRUE)

placenta_seurat <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta/placenta_seurat_v20.RDS")

umap_sct <- Embeddings(placenta_seurat, "umap")
expr_sct <- placenta_seurat@assays$SCT@data["DDX3Y", ]
umap_df <- data.frame(
  UMAP_1 = umap_sct[,1],
  UMAP_2 = umap_sct[,2],
  DDX3Y  = expr_sct,
  stringsAsFactors = FALSE
)
umap_df <- umap_df[order(umap_df$DDX3Y), ]

p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = DDX3Y)) +
  geom_point(size = 0.25) +
  scale_color_gradientn(colors = c("lightgrey", "orange", "purple", "darkblue"), values = c(0, 0.25, 0.5, 1), na.value = "lightgrey") +
  theme_minimal() +
  labs(title = "DDX3Y", x = "UMAP 1", y = "UMAP 2", color = "DDX3Y")

ggsave(filename = file.path(outFolder, "DDX3Y_umap.pdf"), plot = p, device = cairo_pdf, width = 7, height = 6, units = "in", dpi = 600)
ggsave(filename = file.path(outFolder, "DDX3Y_umap.tiff"), plot = p, device = "tiff", width = 7, height = 6, units = "in", dpi = 600)

meta <- placenta_seurat@meta.data
umap_df_split <- data.frame(
  UMAP_1  = umap_sct[,1],
  UMAP_2  = umap_sct[,2],
  DDX3Y   = expr_sct,
  fet_sex = meta[rownames(umap_sct), "fet_sex"],
  stringsAsFactors = FALSE
)
umap_df_split <- umap_df_split[order(umap_df_split$DDX3Y), ]

p_split <- ggplot(umap_df_split, aes(UMAP_1, UMAP_2, color = DDX3Y)) +
  geom_point(size = 0.25) +
  scale_color_gradientn(colors = c("lightgrey", "orange", "purple", "darkblue"), values = c(0, 0.25, 0.5, 1), na.value = "lightgrey") +
  theme_minimal() +
  labs(x = "UMAP 1", y = "UMAP 2", color = "DDX3Y") +
  facet_wrap(~ fet_sex, nrow = 1, labeller = labeller(fet_sex = c("m" = "Male", "f" = "Female")))

ggsave(filename = file.path(outFolder, "DDX3Y_umap_split.pdf"), plot = p_split, device = cairo_pdf, width = 14, height = 6, units = "in", dpi = 600)
ggsave(filename = file.path(outFolder, "DDX3Y_umap_split.tiff"), plot = p_split, device = "tiff", width = 14, height = 6, units = "in", dpi = 600)

umap_df_male   <- subset(umap_df_split, fet_sex == "m")
umap_df_female <- subset(umap_df_split, fet_sex == "f")

make_ddx3y_plot <- function(df, title_txt) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = DDX3Y)) +
    geom_point(size = 0.25) +
    scale_color_gradientn(colors = c("lightgrey", "orange", "purple", "darkblue"), values = c(0, 0.25, 0.5, 1), na.value = "lightgrey") +
    theme_minimal() +
    labs(title = title_txt, x = "UMAP 1", y = "UMAP 2", color = "DDX3Y")
}

p_male   <- make_ddx3y_plot(umap_df_male,   "DDX3Y (Male)")
p_female <- make_ddx3y_plot(umap_df_female, "DDX3Y (Female)")

ggsave(filename = file.path(outFolder, "DDX3Y_umap_male.pdf"),   plot = p_male,   device = cairo_pdf, width = 7, height = 6, units = "in", dpi = 600)
ggsave(filename = file.path(outFolder, "DDX3Y_umap_female.pdf"), plot = p_female, device = cairo_pdf, width = 7, height = 6, units = "in", dpi = 600)

mac_ids <- c("FMAC", "HB", "MAC0", "MAC1", "MAC2", "MAC3", "MAC4", "NEUT") # replace w/ other subclusters
mac <- readRDS("/gscratch/kawaldorflab/jcorn427/placenta/macro_sub_v8.RDS")

mac <- RunPCA(mac, verbose = FALSE)
ElbowPlot(mac, ndims = 50)
mac <- RunUMAP(mac, dims = 1:30, verbose = FALSE)
mac <- FindNeighbors(mac, dims = 1:30, verbose = FALSE)
mac <- FindClusters(mac, resolution = 0.4, verbose = FALSE)

umap_sct_j <- Embeddings(mac, "umap")
expr_sct_j <- mac@assays$SCT@data["DDX3Y", ]
celltype_col <- if ("celltype19" %in% colnames(mac@meta.data)) "celltype19" else if ("celltype18" %in% colnames(mac@meta.data)) "celltype18" else "seurat_clusters"

umap_df_j <- data.frame(
  UMAP_1    = umap_sct_j[,1],
  UMAP_2    = umap_sct_j[,2],
  DDX3Y     = expr_sct_j,
  celltype  = mac@meta.data[rownames(umap_sct_j), celltype_col],
  stringsAsFactors = FALSE
)
