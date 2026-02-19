# =============================================================================
# hdWGCNA gene-centric UMAP visualization
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Gene-centric UMAP of topological overlap with top 10 hub genes per module.
# =============================================================================

library(Seurat)
library(hdWGCNA)
library(ggplot2)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading hdWGCNA object...\n")
mesophyll <- readRDS(snakemake@input[["hdwgcna_obj"]])

cat("Plotting gene-centric UMAP...\n")
p <- PlotModuleUMAP(mesophyll) +
  ggtitle("hdWGCNA Gene-Centric UMAP\n(Topological Overlap with Hub Genes)") +
  theme_minimal()

ggsave(snakemake@output[["gene_umap"]], plot = p,
       width = 8, height = 6, device = "pdf")

cat("Done.\n")
sink(type = "message")
sink(type = "output")
