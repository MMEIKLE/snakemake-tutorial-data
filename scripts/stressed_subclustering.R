# =============================================================================
# Stressed cell subclustering
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Subset the stressed cell type cluster from the full snRNA-seq object,
# re-cluster with Louvain at resolution 0.25, identify subcluster markers,
# and preserve the original UMAP coordinates.
# =============================================================================

library(Seurat)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading snRNA-seq Seurat object...\n")
seurat_obj <- readRDS(snakemake@input[["seurat_obj"]])

# Identify the stressed cell cluster
# The paper describes a "stressed" cell type identified during initial clustering.
# Here we subset cells annotated as stressed (user must set this label during
# annotation, or identify the cluster number corresponding to stressed cells).
cat("Subsetting stressed cells...\n")

# Store original UMAP coordinates before subsetting
original_umap <- Embeddings(seurat_obj, "umap")

# Subset stressed cells — adjust cluster identity as needed
stressed <- subset(seurat_obj, idents = "stressed")

cat(sprintf("  Stressed cells: %d\n", ncol(stressed)))

# Re-cluster the stressed subset at resolution 0.25
cat("Re-clustering stressed cells...\n")
n_pcs <- snakemake@config[["n_pcs"]]

stressed <- FindNeighbors(stressed, dims = 1:n_pcs, verbose = FALSE)
stressed <- FindClusters(
  stressed,
  resolution = 0.25,
  algorithm = 1,  # Louvain
  verbose = FALSE
)

# Preserve original UMAP coordinates
stressed_cells <- colnames(stressed)
stressed[["original_umap"]] <- CreateDimReducObject(
  embeddings = original_umap[stressed_cells, ],
  key = "origUMAP_",
  assay = DefaultAssay(stressed)
)

# Identify subcluster marker genes (Wilcoxon Rank Sum, adj. P < 0.05)
cat("Finding subcluster markers...\n")
subcluster_markers <- FindAllMarkers(
  stressed,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
subcluster_markers <- subcluster_markers[subcluster_markers$p_val_adj < 0.05, ]

stressed@misc$subcluster_markers <- subcluster_markers

cat("Saving stressed subclusters...\n")
saveRDS(stressed, file = snakemake@output[["subclusters"]])

cat(sprintf("  Subclusters: %d\n", length(unique(Idents(stressed)))))
cat(sprintf("  Marker genes: %d\n", nrow(subcluster_markers)))

sink(type = "message")
sink(type = "output")
