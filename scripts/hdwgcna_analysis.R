# =============================================================================
# hdWGCNA co-expression network analysis (v0.3.03)
# Methods: Hu et al. (2025) Front. Plant Sci. 16:1637176
#
# Steps:
#   1. Filter to soybean genes expressed in >=1% of cells
#   2. Subset to healthy and stressed mesophyll cells
#   3. TestSoftPowers for soft-thresholding power selection
#   4. ConstructNetwork: signed, minModuleSize=50,
#      detectCutHeight=0.995, mergeCutHeight=0.2
#   5. Module eigengenes (MEs) and hub genes
#   6. Gene-centric UMAP on topological overlap
#   7. Wilcoxon Rank Sum test for ME comparison between cell states
# =============================================================================

library(Seurat)
library(hdWGCNA)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

cat("Loading snRNA-seq Seurat object...\n")
seurat_obj <- readRDS(snakemake@input[["seurat_obj"]])

# Parameters from config
min_module_size    <- snakemake@config[["hdwgcna_min_module_size"]]     # 50
detect_cut_height  <- snakemake@config[["hdwgcna_detect_cut_height"]]  # 0.995
merge_cut_height   <- snakemake@config[["hdwgcna_merge_cut_height"]]   # 0.2
min_cell_pct       <- snakemake@config[["hdwgcna_min_cell_pct"]]       # 0.01

# --- Step 1-2: Filter genes and subset cells ---
cat("Subsetting to healthy and stressed mesophyll cells...\n")

# Subset to healthy and stressed mesophyll cells
# Adjust cell type labels to match your annotation
mesophyll <- subset(seurat_obj,
  cells = colnames(seurat_obj)[
    seurat_obj$cell_type %in% c("healthy_mesophyll", "stressed_mesophyll")
  ]
)

cat(sprintf("  Cells after subsetting: %d\n", ncol(mesophyll)))

# Set up hdWGCNA
mesophyll <- SetupForWGCNA(
  mesophyll,
  gene_select = "fraction",
  fraction = min_cell_pct,  # genes in >= 1% of cells
  wgcna_name = "hdWGCNA"
)

# Construct metacells for hdWGCNA
mesophyll <- MetacellsByGroups(
  seurat_obj = mesophyll,
  group.by = c("cell_type", "sample"),
  k = 25,
  max_shared = 10,
  ident.group = "cell_type"
)

mesophyll <- NormalizeMetacells(mesophyll)

cat(sprintf("  Genes for network: %d\n",
    length(GetWGCNAGenes(mesophyll, "hdWGCNA"))))

# --- Step 3: Soft power threshold selection ---
cat("Testing soft powers...\n")
mesophyll <- TestSoftPowers(
  mesophyll,
  networkType = "signed"
)

# Select optimal soft power (first power where R^2 > 0.8)
power_table <- GetPowerTable(mesophyll)
soft_power <- min(power_table$Power[power_table$SFT.R.sq > 0.8])
cat(sprintf("  Selected soft power: %d\n", soft_power))

# --- Step 4: Construct network ---
cat("Constructing co-expression network...\n")
mesophyll <- ConstructNetwork(
  mesophyll,
  soft_power = soft_power,
  setDatExpr = FALSE,
  networkType = "signed",
  minModuleSize = min_module_size,
  detectCutHeight = detect_cut_height,
  mergeCutHeight = merge_cut_height,
  tom_name = "hdWGCNA_TOM"
)

cat(sprintf("  Modules identified: %d\n",
    length(unique(GetModules(mesophyll)$module)) - 1))  # exclude grey

# --- Step 5: Module eigengenes and hub genes ---
cat("Computing module eigengenes...\n")
mesophyll <- ModuleEigengenes(mesophyll, group.by.vars = "sample")

# Hub genes: top genes correlated with module eigengene
cat("Identifying hub genes...\n")
mesophyll <- ModuleConnectivity(mesophyll)

modules <- GetModules(mesophyll)
hub_genes <- GetHubGenes(mesophyll, n_hubs = 10)

cat(sprintf("  Hub genes (top 10 per module): %d\n", nrow(hub_genes)))

# --- Step 6: Gene-centric UMAP ---
cat("Running gene-centric UMAP on topological overlap...\n")
mesophyll <- RunModuleUMAP(
  mesophyll,
  n_hubs = 10,
  n_neighbors = 15,
  min_dist = 0.1
)

# --- Step 7: Wilcoxon Rank Sum test on module eigengenes ---
cat("Statistical comparison of module eigengenes between cell states...\n")
me_df <- GetMEs(mesophyll, harmonized = TRUE)
me_df$cell_type <- mesophyll$cell_type

# Test each module eigengene between healthy and stressed
module_names <- colnames(me_df)[grepl("^ME", colnames(me_df))]
wilcox_results <- data.frame(
  module = character(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (mod in module_names) {
  healthy_vals  <- me_df[[mod]][me_df$cell_type == "healthy_mesophyll"]
  stressed_vals <- me_df[[mod]][me_df$cell_type == "stressed_mesophyll"]

  if (length(healthy_vals) > 0 && length(stressed_vals) > 0) {
    wt <- wilcox.test(healthy_vals, stressed_vals, alternative = "two.sided")
    wilcox_results <- rbind(wilcox_results, data.frame(
      module = mod,
      statistic = wt$statistic,
      p_value = wt$p.value
    ))
  }
}

# Adjust p-values
wilcox_results$p_adj <- p.adjust(wilcox_results$p_value, method = "BH")

cat(sprintf("  Significant modules (padj < 0.05): %d\n",
    sum(wilcox_results$p_adj < 0.05)))

# --- Save outputs ---
cat("Saving results...\n")
saveRDS(mesophyll, file = snakemake@output[["hdwgcna_obj"]])

# Export module eigengenes
me_export <- GetMEs(mesophyll, harmonized = TRUE)
me_export$cell <- rownames(me_export)
write.csv(me_export, file = snakemake@output[["eigengenes"]], row.names = FALSE)

# Export hub genes
write.csv(hub_genes, file = snakemake@output[["hub_genes"]], row.names = FALSE)

cat("hdWGCNA analysis complete.\n")
sink(type = "message")
sink(type = "output")
