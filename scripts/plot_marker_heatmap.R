# =============================================================================
# Marker gene heatmaps using pheatmap (v1.0.12)
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
# =============================================================================

library(Seurat)
library(pheatmap)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading data...\n")
spatial_markers <- read.csv(snakemake@input[["spatial_markers"]])
snrna_markers   <- read.csv(snakemake@input[["snrna_markers"]])
spatial_obj     <- readRDS(snakemake@input[["spatial_obj"]])
snrna_obj       <- readRDS(snakemake@input[["snrna_obj"]])

# Select top 5 markers per cluster for heatmap
get_top_markers <- function(markers_df, n = 5) {
  markers_df <- markers_df[order(markers_df$avg_log2FC, decreasing = TRUE), ]
  top <- do.call(rbind, lapply(split(markers_df, markers_df$cluster), function(x) {
    head(x, n)
  }))
  unique(top$gene)
}

cat("Generating heatmap...\n")

# Combine top markers from both datasets
spatial_top <- get_top_markers(spatial_markers)
snrna_top   <- get_top_markers(snrna_markers)

# Use snRNA-seq object for the heatmap (more cell types)
DefaultAssay(snrna_obj) <- "SCT"

# Average expression per cluster
avg_expr <- AverageExpression(snrna_obj, features = snrna_top,
                               assays = "SCT", return.seurat = FALSE)
mat <- as.matrix(avg_expr$SCT)

# Scale rows for visualization
mat_scaled <- t(scale(t(mat)))

pdf(snakemake@output[["heatmap"]], width = 10, height = 12)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Top Cluster Marker Genes"
)
dev.off()

cat("Done.\n")
sink(type = "message")
sink(type = "output")
