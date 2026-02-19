# =============================================================================
# Seurat snRNA-seq analysis
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Steps:
#   1. Load Cell Ranger outputs for all 6 snRNA-seq samples
#   2. QC filter: 300 < nCount_RNA < 10000; 250 < nFeature_RNA < 6000;
#      percent.mt + percent.chloro < 10%
#   3. DoubletFinder (v2.0.4) to identify and remove doublets
#   4. SCTransform normalization per sample
#   5. CCA integration across all samples
#   6. PCA (30 PCs), Louvain clustering (resolution=0.25), UMAP (30 PCs)
# =============================================================================

library(Seurat)
library(hdf5r)
library(DoubletFinder)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# QC thresholds from the paper
min_umi    <- snakemake@config[["snrna_min_umi"]]     # 300
max_umi    <- snakemake@config[["snrna_max_umi"]]     # 10000
min_genes  <- snakemake@config[["snrna_min_genes"]]   # 250
max_genes  <- snakemake@config[["snrna_max_genes"]]   # 6000
max_mt_pct <- snakemake@config[["snrna_max_mito_pct"]] # 10

sample_paths <- snakemake@input[["matrices"]]
sample_names <- snakemake@config[["snrna_samples"]]

cat("Loading snRNA-seq samples...\n")

seurat_list <- lapply(seq_along(sample_paths), function(i) {
  mat <- Read10X_h5(sample_paths[i])
  obj <- CreateSeuratObject(counts = mat, project = sample_names[i])
  obj$sample <- sample_names[i]
  obj
})
names(seurat_list) <- sample_names

cat("QC filtering...\n")

# Calculate mitochondrial and chloroplast percentages
# Soybean mitochondrial genes typically on chromosome Mt, chloroplast on Pt
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-|^GlymaM")
  obj[["percent.chloro"]] <- PercentageFeatureSet(obj, pattern = "^Pt-|^GlymaPt")
  obj[["percent.organelle"]] <- obj$percent.mt + obj$percent.chloro

  # Apply QC filters
  obj <- subset(obj,
    subset = nCount_RNA > min_umi &
             nCount_RNA < max_umi &
             nFeature_RNA > min_genes &
             nFeature_RNA < max_genes &
             percent.organelle < max_mt_pct
  )
  obj
})

cat("DoubletFinder doublet removal...\n")

# DoubletFinder (v2.0.4) on each sample
seurat_list <- lapply(seurat_list, function(obj) {
  # Pre-process for DoubletFinder
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                              verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.25, verbose = FALSE)

  # Estimate doublet rate (~0.8% per 1000 cells loaded)
  n_cells <- ncol(obj)
  doublet_rate <- n_cells * 0.8e-3 / 1000 * n_cells
  doublet_rate <- min(doublet_rate / n_cells, 0.1)  # cap at 10%

  # Parameter sweep for pK
  sweep_res <- paramSweep(obj, PCs = 1:30, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  pk_opt <- find.pK(sweep_stats)
  pk <- as.numeric(as.character(pk_opt$pK[which.max(pk_opt$BCmetric)]))

  # Estimate homotypic doublet proportion
  annotations <- obj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  nExp <- round(doublet_rate * n_cells)
  nExp_adj <- round(nExp * (1 - homotypic_prop))

  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:30, pN = 0.25, pK = pk,
                       nExp = nExp_adj, sct = FALSE)

  # Remove doublets
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
  obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[df_col]] == "Singlet"])
  obj
})

cat("SCTransform normalization per sample...\n")

seurat_list <- lapply(seurat_list, function(obj) {
  SCTransform(obj, verbose = FALSE)
})

cat("CCA integration...\n")

features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)
seurat_list <- PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = features
)

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

n_pcs <- snakemake@config[["n_pcs"]]
resolution <- snakemake@config[["snrna_cluster_resolution"]]

integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:n_pcs, verbose = FALSE)
integrated <- FindClusters(
  integrated,
  resolution = resolution,
  algorithm = 1,  # Louvain
  verbose = FALSE
)
integrated <- RunUMAP(integrated, dims = 1:n_pcs, verbose = FALSE)

cat("Saving Seurat object...\n")
saveRDS(integrated, file = snakemake@output[["seurat_obj"]])

cat("snRNA-seq analysis complete.\n")
cat(sprintf("  Nuclei: %d\n", ncol(integrated)))
cat(sprintf("  Genes: %d\n", nrow(integrated)))
cat(sprintf("  Clusters: %d\n", length(unique(Idents(integrated)))))

sink(type = "message")
sink(type = "output")
