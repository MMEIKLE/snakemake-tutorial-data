# =============================================================================
# Pseudo-bulk analysis with DESeq2 (v1.42.1)
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Aggregate raw counts per gene using Seurat AggregateExpression,
# then normalize with DESeq2 default size factor estimation.
# =============================================================================

library(Seurat)
library(DESeq2)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(snakemake@input[["seurat_obj"]])

# Aggregate raw counts per gene, grouped by sample
cat("Aggregating expression by sample...\n")
DefaultAssay(seurat_obj) <- "RNA"

pseudo_bulk <- AggregateExpression(
  seurat_obj,
  group.by = "sample",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

count_matrix <- as.matrix(pseudo_bulk$RNA)

cat(sprintf("  Genes: %d\n", nrow(count_matrix)))
cat(sprintf("  Samples: %d\n", ncol(count_matrix)))

# DESeq2 normalization
cat("Running DESeq2 normalization...\n")
col_data <- data.frame(
  sample = colnames(count_matrix),
  row.names = colnames(count_matrix)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),  # ensure integer counts
  colData = col_data,
  design = ~ 1  # normalization only, no differential testing
)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

cat(sprintf("  Size factors: %s\n",
    paste(round(sizeFactors(dds), 3), collapse = ", ")))

write.csv(normalized_counts, file = snakemake@output[["normalized"]])

cat("Done.\n")
sink(type = "message")
sink(type = "output")
