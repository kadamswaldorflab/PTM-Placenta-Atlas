library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

outFolder <- "./HumanNHP_mmulAnalysis/Final_Pipeline/GOanalysis/"
dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)

placenta_seurat <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta/placenta_seurat_v20.RDS")
Idents(placenta_seurat) <- placenta_seurat$celltype19
clusters <- levels(Idents(placenta_seurat))

up_genes <- list()
down_genes <- list()

for (cl in clusters) {
  m_out <- FindMarkers(
    object = placenta_seurat,
    ident.1 = cl,
    ident.2 = NULL,
    only.pos = FALSE,
    logfc.threshold = 0
  )
  top_up <- m_out %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 3) %>%
    rownames()
  top_down <- m_out %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 3) %>%
    rownames()
  up_genes[[cl]] <- top_up
  down_genes[[cl]] <- top_down
}

cc_up <- compareCluster(
  geneCluster = up_genes,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  keyType = "SYMBOL",
  readable = TRUE
)

cc_down <- compareCluster(
  geneCluster = down_genes,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  keyType = "SYMBOL",
  readable = TRUE
)

saveRDS(cc_up, file.path(outFolder, "cc_up.rds"))
saveRDS(cc_down, file.path(outFolder, "cc_down.rds"))

subset_cc <- function(cc_obj, keep_clusters) {
  cc2 <- cc_obj
  cc2@compareClusterResult <- cc_obj@compareClusterResult %>% filter(Cluster %in% keep_clusters)
  cc2
}

group_list <- list(
  CTB = c("schCTB", "CTB", "iCTB", "STB", "PC", "FIB", "EC1", "EC2"),
  Immune = c("uNKP", "uNK1", "uNK2", "NKT", "CD4+Treg", "TRM CD8+T", "IgM+B", "IgG+B", "Mast"),
  Macs = c("FMAC", "HB", "MAC0", "MAC1", "MAC2", "MAC3", "MAC4", "NEUT"),
  EB = c("EB1", "EB2"),
  DSC = c("DMSC", "DSC0", "DSC1", "DSC2", "DSC3", "DSC4")
)

for (grp_name in names(group_list)) {
  cls <- group_list[[grp_name]]
  cc_sub_up <- subset_cc(cc_up, cls)
  p_up <- dotplot(cc_sub_up, showCategory = 3) + ggtitle(paste0("Up-regulated GO (BP): ", grp_name))
  ggsave(filename = file.path(outFolder, paste0("BP_dotplot_up_", grp_name, ".pdf")), plot = p_up, width = 10, height = 15)
  cc_sub_dn <- subset_cc(cc_down, cls)
  p_dn <- dotplot(cc_sub_dn, showCategory = 3) + ggtitle(paste0("Down-regulated GO (BP): ", grp_name))
  ggsave(filename = file.path(outFolder, paste0("BP_dotplot_down_", grp_name, ".pdf")), plot = p_dn, width = 10, height = 15)
}
