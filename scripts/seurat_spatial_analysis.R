# =============================================================================
# Seurat spatial transcriptomics analysis (Visium)
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Steps:
#   1. Load Space Ranger outputs for all 12 Visium samples
#   2. SCTransform normalization per sample (Hafemeister & Satija, 2019)
#   3. CCA integration across all capture areas (Seurat v5.1.0)
#   4. PCA (30 PCs), Louvain clustering (resolution=0.5), UMAP (30 PCs)
# =============================================================================

library(Seurat)
library(hdf5r)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading Visium samples...\n")

# Load all Visium samples
sample_paths <- snakemake@input[["matrices"]]
sample_names <- snakemake@config[["visium_samples"]]

seurat_list <- lapply(seq_along(sample_paths), function(i) {
  sample_dir <- dirname(dirname(sample_paths[i]))  # up to {sample}/outs
  obj <- Load10X_Spatial(
    data.dir = file.path(sample_dir, "outs"),
    filename = "filtered_feature_bc_matrix.h5"
  )
  obj$sample <- sample_names[i]
  obj
})
names(seurat_list) <- sample_names

cat("SCTransform normalization per sample...\n")

# SCTransform normalization per sample
seurat_list <- lapply(seurat_list, function(obj) {
  SCTransform(obj, assay = "Spatial", verbose = FALSE)
})

cat("CCA integration across all capture areas...\n")

# Select integration features and prepare for CCA
features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)
seurat_list <- PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = features
)

# CCA integration
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "cca"
)
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

cat("PCA, clustering, and UMAP...\n")

# PCA
n_pcs <- snakemake@config[["n_pcs"]]
integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)

# Louvain clustering at resolution 0.5
resolution <- snakemake@config[["spatial_cluster_resolution"]]
integrated <- FindNeighbors(integrated, dims = 1:n_pcs, verbose = FALSE)
integrated <- FindClusters(
  integrated,
  resolution = resolution,
  algorithm = 1,  # Louvain
  verbose = FALSE
)

# UMAP
integrated <- RunUMAP(integrated, dims = 1:n_pcs, verbose = FALSE)

cat("Saving Seurat object...\n")
saveRDS(integrated, file = snakemake@output[["seurat_obj"]])

cat("Spatial analysis complete.\n")
cat(sprintf("  Spots: %d\n", ncol(integrated)))
cat(sprintf("  Genes: %d\n", nrow(integrated)))
cat(sprintf("  Clusters: %d\n", length(unique(Idents(integrated)))))

sink(type = "message")
sink(type = "output")
