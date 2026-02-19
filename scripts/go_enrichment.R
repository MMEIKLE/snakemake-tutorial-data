# =============================================================================
# Gene Ontology enrichment with fgsea (v1.28.0)
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Genes ranked by avg_log2FC from FindAllMarkers or FindMarkers output.
# GO gene sets provided in GMT format.
# =============================================================================

library(fgsea)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading marker genes...\n")
markers <- read.csv(snakemake@input[["markers"]])

cat("Loading GO gene sets...\n")
go_pathways <- gmtPathways(snakemake@input[["go_genesets"]])

# Build ranked gene list by avg_log2FC
# If multiple entries per gene (from FindAllMarkers across clusters),
# use the maximum log2FC
gene_ranks <- tapply(markers$avg_log2FC, markers$gene, max)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat(sprintf("  Ranked genes: %d\n", length(gene_ranks)))
cat(sprintf("  GO pathways: %d\n", length(go_pathways)))

# Run fgsea
cat("Running fgsea...\n")
fgsea_res <- fgsea(
  pathways = go_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500
)

# Filter significant results
fgsea_sig <- fgsea_res[fgsea_res$padj < 0.05, ]
fgsea_sig <- fgsea_sig[order(fgsea_sig$padj), ]

cat(sprintf("  Significant GO terms (padj < 0.05): %d\n", nrow(fgsea_sig)))

# Convert leadingEdge list column to semicolon-separated string for CSV
fgsea_sig$leadingEdge <- sapply(fgsea_sig$leadingEdge, paste, collapse = ";")

write.csv(
  as.data.frame(fgsea_sig),
  file = snakemake@output[["enrichment"]],
  row.names = FALSE
)

cat("Done.\n")
sink(type = "message")
sink(type = "output")
