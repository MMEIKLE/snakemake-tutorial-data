# =============================================================================
# Venn diagrams and UpSet plots
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# eulerr (v7.0.2) + SuperExactTest (v1.1.0) for Venn diagrams
# UpSetR (v1.4.0) for UpSet plots
# =============================================================================

library(eulerr)
library(SuperExactTest)
library(UpSetR)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading marker data...\n")
spatial_markers <- read.csv(snakemake@input[["spatial_markers"]])
snrna_markers   <- read.csv(snakemake@input[["snrna_markers"]])

# Build gene sets per cluster
spatial_sets <- split(spatial_markers$gene, spatial_markers$cluster)
snrna_sets   <- split(snrna_markers$gene, snrna_markers$cluster)

# --- Venn diagram: overlap between spatial and snRNA-seq marker genes ---
cat("Creating Venn diagram...\n")

all_spatial_genes <- unique(spatial_markers$gene)
all_snrna_genes   <- unique(snrna_markers$gene)

venn_input <- list(
  Spatial = all_spatial_genes,
  snRNA = all_snrna_genes
)

# SuperExactTest for statistical significance of overlap
total_genes <- length(union(all_spatial_genes, all_snrna_genes))
set_test <- supertest(venn_input, n = total_genes)

pdf(snakemake@output[["venn"]], width = 6, height = 6)

# eulerr proportional Venn
fit <- euler(venn_input)
plot(fit,
  fills = list(fill = c("#E41A1C", "#377EB8"), alpha = 0.5),
  labels = list(fontsize = 12),
  quantities = list(fontsize = 10),
  main = sprintf("Marker Gene Overlap\n(SuperExactTest p = %.2e)",
                  set_test$P.value)
)

dev.off()

# --- UpSet plot: overlap across individual clusters ---
cat("Creating UpSet plot...\n")

# Combine spatial and snRNA cluster gene sets
all_sets <- c(
  setNames(spatial_sets, paste0("Sp_", names(spatial_sets))),
  setNames(snrna_sets, paste0("sn_", names(snrna_sets)))
)

# Limit to top clusters by marker count for readability
set_sizes <- sapply(all_sets, length)
top_sets <- names(sort(set_sizes, decreasing = TRUE))[1:min(10, length(all_sets))]
all_sets <- all_sets[top_sets]

pdf(snakemake@output[["upset"]], width = 12, height = 8)
upset(
  fromList(all_sets),
  order.by = "freq",
  nsets = length(all_sets),
  text.scale = 1.2,
  mainbar.y.label = "Shared Marker Genes",
  sets.x.label = "Markers per Cluster"
)
dev.off()

cat("Done.\n")
sink(type = "message")
sink(type = "output")
