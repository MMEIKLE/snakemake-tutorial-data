# =============================================================================
# Find cluster marker genes — spatial data
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Seurat FindAllMarkers with Wilcoxon Rank Sum test, adjusted P < 0.05.
# =============================================================================

library(Seurat)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(snakemake@input[["seurat_obj"]])

# Ensure we use the SCT assay for marker detection
DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- PrepSCTFindMarkers(seurat_obj)

cat("Running FindAllMarkers (Wilcoxon Rank Sum)...\n")
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Filter to adjusted P < 0.05
markers <- markers[markers$p_val_adj < 0.05, ]

cat(sprintf("  Significant markers: %d\n", nrow(markers)))
cat(sprintf("  Clusters with markers: %d\n", length(unique(markers$cluster))))

write.csv(markers, file = snakemake@output[["markers"]], row.names = FALSE)

cat("Done.\n")
sink(type = "message")
sink(type = "output")
