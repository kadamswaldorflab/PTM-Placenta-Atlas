library(Seurat)
library(sccomp)
library(tidyverse)

seurat_obj <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/placenta/placenta_seurat_v20.RDS")

seurat_obj$Tissue <- sub(".*([A-Za-z])$", "\\1", seurat_obj$sample_id)
seurat_obj$Subgroup <- ifelse(grepl("^CTRL", seurat_obj$sample_id), "CTRL",
                       ifelse(grepl("^SAL", seurat_obj$sample_id), "SAL", NA))

metadata_fetsex <- seurat_obj@meta.data %>%
  as_tibble(rownames = "cell_barcode") %>%
  transmute(
    sample_id,
    fetsex = recode(as.character(fet_sex), "f" = "F", "m" = "M"),
    cell_group = celltype19
  )

estimate_fetsex <- sccomp_estimate(
  .data = metadata_fetsex,
  formula_composition = ~ fetsex,
  formula_variability = ~ 1,
  .sample = sample_id,
  .cell_group = cell_group,
  cores = 1
)

estimate_fetsex_tested <- estimate_fetsex %>% sccomp_test()

sex_sccomp2 <- sccomp_boxplot(
  .data = estimate_fetsex_tested,
  factor = "fetsex",
  significance_threshold = 0.1
) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.minor = element_blank()
  )

print(sex_sccomp2)

metadata_subgroup <- seurat_obj@meta.data %>%
  as_tibble(rownames = "cell_barcode") %>%
  select(
    sample_id,
    Subgroup,
    cell_group = celltype19
  )

estimate_subgroup <- sccomp_estimate(
  .data = metadata_subgroup,
  formula_composition = ~ Subgroup,
  formula_variability = ~ 1,
  .sample = sample_id,
  .cell_group = cell_group,
  .abundance = NULL,
  cores = 1
)

estimate_subgroup_tested <- estimate_subgroup %>% sccomp_test()

subgroup_sccomp <- sccomp_boxplot(
  .data = estimate_subgroup_tested,
  factor = "Subgroup",
  significance_threshold = 0.05
) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.minor = element_blank()
  )

print(subgroup_sccomp)

metadata_tissue <- seurat_obj@meta.data %>%
  as_tibble(rownames = "cell_barcode") %>%
  select(
    sample_id,
    Tissue,
    cell_group = celltype19
  ) %>%
  mutate(
    Tissue = recode(Tissue, "c" = "DISK", "m" = "CAM", "d" = "DEC")
  )

estimate_tissue <- sccomp_estimate(
  .data = metadata_tissue,
  formula_composition = ~ Tissue,
  formula_variability = ~ 1,
  .sample = sample_id,
  .cell_group = cell_group,
  .abundance = NULL,
  cores = 1
)

estimate_tissue_tested <- estimate_tissue %>% sccomp_test()

tissue_sccomp <- sccomp_boxplot(
  .data = estimate_tissue_tested,
  factor = "Tissue",
  significance_threshold = 0.05
) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.minor = element_blank()
  )

print(tissue_sccomp)
