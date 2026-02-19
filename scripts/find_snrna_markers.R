# =============================================================================
# Find cluster markers and DEGs — snRNA-seq data
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# 1. FindAllMarkers: cluster markers (Wilcoxon Rank Sum, adj. P < 0.05)
# 2. FindMarkers: healthy vs infected DEGs (Wilcoxon Rank Sum, adj. P < 0.05)
# 3. AddModuleScore: 15 ASR-induced genes (Cabre et al. 2021; Morales 2013)
# =============================================================================

library(Seurat)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(snakemake@input[["seurat_obj"]])

DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- PrepSCTFindMarkers(seurat_obj)

# --- FindAllMarkers: cluster markers ---
cat("Running FindAllMarkers (Wilcoxon Rank Sum)...\n")
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
all_markers <- all_markers[all_markers$p_val_adj < 0.05, ]

cat(sprintf("  Cluster markers: %d\n", nrow(all_markers)))
write.csv(all_markers, file = snakemake@output[["markers"]], row.names = FALSE)

# --- FindMarkers: healthy vs infected DEGs ---
# The paper compares healthy and infected clusters. Adjust ident.1/ident.2
# to match your cluster annotations.
cat("Running FindMarkers: healthy vs infected...\n")
degs <- FindMarkers(
  seurat_obj,
  ident.1 = "infected",
  ident.2 = "healthy",
  test.use = "wilcox",
  logfc.threshold = 0,  # return all genes for fgsea ranking
  min.pct = 0.1
)
degs$gene <- rownames(degs)
degs <- degs[degs$p_val_adj < 0.05, ]

cat(sprintf("  DEGs (healthy vs infected): %d\n", nrow(degs)))
write.csv(degs, file = snakemake@output[["degs"]], row.names = FALSE)

# --- AddModuleScore: ASR-induced gene set ---
cat("Adding ASR-induced gene module score...\n")
asr_genes <- readLines(snakemake@input[["asr_genes"]])
asr_genes <- asr_genes[asr_genes != ""]

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(ASR_induced = asr_genes),
  name = "ASR_module_score"
)

# Save updated object with module scores
saveRDS(seurat_obj, file = snakemake@input[["seurat_obj"]])

cat("Done.\n")
sink(type = "message")
sink(type = "output")
