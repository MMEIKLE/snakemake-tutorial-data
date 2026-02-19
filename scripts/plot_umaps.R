# =============================================================================
# UMAP visualizations for spatial and snRNA-seq data
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
# =============================================================================

library(Seurat)
library(ggplot2)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# --- Spatial UMAP ---
cat("Plotting spatial UMAP...\n")
spatial_obj <- readRDS(snakemake@input[["spatial_obj"]])

p_spatial <- DimPlot(spatial_obj, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, pt.size = 0.5) +
  ggtitle("Spatial Transcriptomics - UMAP") +
  theme_minimal()

ggsave(snakemake@output[["spatial_umap"]], plot = p_spatial,
       width = 8, height = 6, device = "pdf")

# --- snRNA-seq UMAP ---
cat("Plotting snRNA-seq UMAP...\n")
snrna_obj <- readRDS(snakemake@input[["snrna_obj"]])

p_snrna <- DimPlot(snrna_obj, reduction = "umap", group.by = "seurat_clusters",
                   label = TRUE, pt.size = 0.3) +
  ggtitle("snRNA-seq - UMAP") +
  theme_minimal()

ggsave(snakemake@output[["snrna_umap"]], plot = p_snrna,
       width = 8, height = 6, device = "pdf")

cat("Done.\n")
sink(type = "message")
sink(type = "output")
