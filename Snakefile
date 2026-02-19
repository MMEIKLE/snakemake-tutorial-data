# =============================================================================
# Snakemake pipeline: Spatial and single-cell transcriptomics of soybean
# defense response to Phakopsora pachyrhizi (Asian Soybean Rust)
#
# Reproduces the computational methods from:
#   Hu et al. (2025) "Spatial and single-cell transcriptomics capture two
#   distinct cell states in soybean defense response to Phakopsora pachyrhizi
#   infection." Front. Plant Sci. 16:1637176.
#   doi: 10.3389/fpls.2025.1637176
#   Data: PRJNA1291134
#
# Pipeline overview:
#   1. Build combined soybean + ASR reference genome
#   2. Space Ranger (v2.0.1) for Visium spatial transcriptomics
#   3. Cell Ranger (v8.0.1) for snRNA-seq
#   4. QC filtering, doublet removal (DoubletFinder v2.0.4)
#   5. SCTransform normalization + CCA integration (Seurat v5.1.0)
#   6. PCA, Louvain clustering, UMAP
#   7. Stressed cell subclustering
#   8. Marker gene identification (Wilcoxon Rank Sum, adj. P < 0.05)
#   9. GO enrichment (fgsea v1.28.0)
#  10. Pseudo-bulk DESeq2 (v1.42.1) normalization
#  11. hdWGCNA co-expression network analysis (v0.3.03)
#  12. Visualizations (pheatmap, eulerr, UpSetR, SuperExactTest)
# =============================================================================

import os

configfile: "config.yaml"


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # Space Ranger outputs
        expand(
            "results/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
            sample=config["visium_samples"],
        ),
        # Cell Ranger outputs
        expand(
            "results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
            sample=config["snrna_samples"],
        ),
        # Seurat analysis outputs
        "results/spatial/seurat_spatial_clustered.rds",
        "results/snrna/seurat_snrna_clustered.rds",
        "results/snrna/stressed_subclusters.rds",
        # Marker genes and DE
        "results/spatial/spatial_cluster_markers.csv",
        "results/snrna/snrna_cluster_markers.csv",
        "results/snrna/healthy_vs_infected_DEGs.csv",
        # GO enrichment
        "results/spatial/spatial_GO_enrichment.csv",
        "results/snrna/snrna_GO_enrichment.csv",
        # Pseudo-bulk
        "results/spatial/pseudobulk_deseq2_normalized.csv",
        "results/snrna/pseudobulk_deseq2_normalized.csv",
        # hdWGCNA
        "results/hdwgcna/hdwgcna_modules.rds",
        "results/hdwgcna/module_eigengenes.csv",
        "results/hdwgcna/hub_genes.csv",
        # Visualizations
        "results/figures/spatial_umap.pdf",
        "results/figures/snrna_umap.pdf",
        "results/figures/marker_heatmap.pdf",
        "results/figures/venn_diagram.pdf",
        "results/figures/upset_plot.pdf",
        "results/figures/hdwgcna_gene_umap.pdf",


# ===========================================================================
# STEP 1: Build combined reference genome (soybean Wm82v4 + ASR Phapa1)
# ===========================================================================
rule build_combined_reference:
    """
    Combine Glycine max (Wm82v4) and Phakopsora pachyrhizi (Phapa1,
    GCA_025201825.1) genomes and gene annotations into a single Cell Ranger-
    compatible reference using cellranger mkref.
    """
    input:
        soy_fasta=config["soybean_genome_fasta"],
        soy_gtf=config["soybean_gtf"],
        asr_fasta=config["asr_genome_fasta"],
        asr_gtf=config["asr_gtf"],
    output:
        directory("reference/combined_soy_asr"),
    log:
        "logs/build_combined_reference.log",
    params:
        mem_gb=config.get("mkref_mem_gb", 64),
    threads: 8
    shell:
        """
        # Concatenate genomes and annotations
        cat {input.soy_fasta} {input.asr_fasta} > reference/combined_genome.fa
        cat {input.soy_gtf} {input.asr_gtf} > reference/combined_genes.gtf

        # Build Cell Ranger reference
        cellranger mkref \
            --genome=combined_soy_asr \
            --fasta=reference/combined_genome.fa \
            --genes=reference/combined_genes.gtf \
            --nthreads={threads} \
            --memgb={params.mem_gb} \
            2>&1 | tee {log}

        mv combined_soy_asr reference/combined_soy_asr
        rm -f reference/combined_genome.fa reference/combined_genes.gtf
        """


# ===========================================================================
# STEP 2: Space Ranger (v2.0.1) — Visium spatial transcriptomics
# ===========================================================================
rule spaceranger_count:
    """
    Process Visium spatial transcriptomics FASTQ files with Space Ranger.
    Aligns reads to the combined soybean + ASR reference and generates
    per-spot gene expression matrices.
    """
    input:
        fastqs=lambda wc: config["visium_fastq_dirs"][wc.sample],
        image=lambda wc: config["visium_images"][wc.sample],
        reference="reference/combined_soy_asr",
    output:
        matrix="results/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        spatial="results/spaceranger/{sample}/outs/spatial/tissue_positions.csv",
    log:
        "logs/spaceranger/{sample}.log",
    params:
        slide=lambda wc: config["visium_slides"][wc.sample],
        area=lambda wc: config["visium_areas"][wc.sample],
    threads: 16
    resources:
        mem_mb=64000,
    shell:
        """
        spaceranger count \
            --id={wildcards.sample} \
            --transcriptome={input.reference} \
            --fastqs={input.fastqs} \
            --image={input.image} \
            --slide={params.slide} \
            --area={params.area} \
            --localcores={threads} \
            --localmem={resources.mem_mb} \
            2>&1 | tee {log}

        # Move outputs to results directory
        mv {wildcards.sample}/outs results/spaceranger/{wildcards.sample}/outs
        rm -rf {wildcards.sample}
        """


# ===========================================================================
# STEP 3: Cell Ranger (v8.0.1) — snRNA-seq
# ===========================================================================
rule cellranger_count:
    """
    Process 10x Chromium snRNA-seq FASTQ files with Cell Ranger.
    Aligns reads to the combined soybean + ASR reference and generates
    per-nucleus gene expression matrices.
    """
    input:
        fastqs=lambda wc: config["snrna_fastq_dirs"][wc.sample],
        reference="reference/combined_soy_asr",
    output:
        matrix="results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
    log:
        "logs/cellranger/{sample}.log",
    threads: 16
    resources:
        mem_mb=64000,
    shell:
        """
        cellranger count \
            --id={wildcards.sample} \
            --transcriptome={input.reference} \
            --fastqs={input.fastqs} \
            --include-introns=true \
            --localcores={threads} \
            --localmem={resources.mem_mb} \
            2>&1 | tee {log}

        mv {wildcards.sample}/outs results/cellranger/{wildcards.sample}/outs
        rm -rf {wildcards.sample}
        """


# ===========================================================================
# STEP 4-6: Spatial Seurat analysis
#   - QC filtering (tissue spots from Loupe Browser alignment)
#   - SCTransform normalization per sample
#   - CCA integration across all capture areas
#   - PCA (30 PCs), Louvain clustering (resolution 0.5), UMAP
# ===========================================================================
rule seurat_spatial_analysis:
    """
    Seurat v5.1.0 workflow for Visium spatial data:
      1. Load Space Ranger outputs for all 12 samples
      2. Filter to tissue-covered spots (pre-selected via Loupe Browser)
      3. SCTransform normalization per sample
      4. CCA integration across all capture areas
      5. PCA (30 PCs), Louvain clustering (resolution=0.5), UMAP (30 PCs)
    """
    input:
        matrices=expand(
            "results/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
            sample=config["visium_samples"],
        ),
    output:
        seurat_obj="results/spatial/seurat_spatial_clustered.rds",
    log:
        "logs/seurat_spatial_analysis.log",
    threads: 4
    resources:
        mem_mb=32000,
    script:
        "scripts/seurat_spatial_analysis.R"


# ===========================================================================
# STEP 4-6: snRNA-seq Seurat analysis
#   - QC: 300 < UMI < 10000, 250 < genes < 6000, <10% mito/chloro
#   - DoubletFinder (v2.0.4) doublet removal
#   - SCTransform normalization per sample
#   - CCA integration across all 6 samples
#   - PCA (30 PCs), Louvain clustering (resolution 0.25), UMAP
# ===========================================================================
rule seurat_snrna_analysis:
    """
    Seurat v5.1.0 workflow for snRNA-seq data:
      1. Load Cell Ranger outputs for all 6 samples
      2. QC filter: 300 < nCount_RNA < 10000; 250 < nFeature_RNA < 6000;
         percent.mt + percent.chloro < 10%
      3. DoubletFinder (v2.0.4) to identify and remove doublets
      4. SCTransform normalization per sample
      5. CCA integration across all samples
      6. PCA (30 PCs), Louvain clustering (resolution=0.25), UMAP (30 PCs)
    """
    input:
        matrices=expand(
            "results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
            sample=config["snrna_samples"],
        ),
    output:
        seurat_obj="results/snrna/seurat_snrna_clustered.rds",
    log:
        "logs/seurat_snrna_analysis.log",
    threads: 4
    resources:
        mem_mb=32000,
    script:
        "scripts/seurat_snrna_analysis.R"


# ===========================================================================
# STEP 7: Stressed cell subclustering
# ===========================================================================
rule stressed_cell_subclustering:
    """
    Subset the stressed cell cluster from the full snRNA-seq Seurat object
    and re-cluster with Louvain at resolution 0.25. Marker genes for
    subclusters are identified and compared to the full dataset clusters.
    Original UMAP coordinates are preserved.
    """
    input:
        seurat_obj="results/snrna/seurat_snrna_clustered.rds",
    output:
        subclusters="results/snrna/stressed_subclusters.rds",
    log:
        "logs/stressed_subclustering.log",
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "scripts/stressed_subclustering.R"


# ===========================================================================
# STEP 8: Marker gene identification
#   - FindAllMarkers: cluster markers (Wilcoxon Rank Sum, adj. P < 0.05)
#   - FindMarkers: healthy vs infected DEGs (Wilcoxon Rank Sum, adj. P < 0.05)
#   - AddModuleScore: ASR-induced gene set (Cabre et al. 2021, Morales 2013)
# ===========================================================================
rule find_spatial_markers:
    """
    Identify cluster marker genes in spatial data using Seurat FindAllMarkers
    with Wilcoxon Rank Sum test (adjusted P < 0.05).
    """
    input:
        seurat_obj="results/spatial/seurat_spatial_clustered.rds",
    output:
        markers="results/spatial/spatial_cluster_markers.csv",
    log:
        "logs/find_spatial_markers.log",
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "scripts/find_markers.R"


rule find_snrna_markers:
    """
    Identify cluster marker genes in snRNA-seq data using Seurat
    FindAllMarkers with Wilcoxon Rank Sum test (adjusted P < 0.05).
    Also runs FindMarkers for healthy vs infected DEGs and AddModuleScore
    for the 15 ASR-induced genes from Cabre et al. (2021) and Morales et al.
    (2013).
    """
    input:
        seurat_obj="results/snrna/seurat_snrna_clustered.rds",
        asr_genes=config["asr_induced_genes"],
    output:
        markers="results/snrna/snrna_cluster_markers.csv",
        degs="results/snrna/healthy_vs_infected_DEGs.csv",
    log:
        "logs/find_snrna_markers.log",
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "scripts/find_snrna_markers.R"


# ===========================================================================
# STEP 9: GO enrichment with fgsea (v1.28.0)
#   Ranked by avg_log2FC from FindAllMarkers / FindMarkers
# ===========================================================================
rule go_enrichment_spatial:
    """
    Gene Ontology enrichment analysis on spatial cluster markers using
    fgsea (v1.28.0). Genes ranked by avg_log2FC.
    """
    input:
        markers="results/spatial/spatial_cluster_markers.csv",
        go_genesets=config["go_genesets"],
    output:
        enrichment="results/spatial/spatial_GO_enrichment.csv",
    log:
        "logs/go_enrichment_spatial.log",
    threads: 2
    script:
        "scripts/go_enrichment.R"


rule go_enrichment_snrna:
    """
    Gene Ontology enrichment analysis on snRNA-seq DEGs using fgsea (v1.28.0).
    Genes ranked by avg_log2FC.
    """
    input:
        markers="results/snrna/healthy_vs_infected_DEGs.csv",
        go_genesets=config["go_genesets"],
    output:
        enrichment="results/snrna/snrna_GO_enrichment.csv",
    log:
        "logs/go_enrichment_snrna.log",
    threads: 2
    script:
        "scripts/go_enrichment.R"


# ===========================================================================
# STEP 10: Pseudo-bulk analysis with DESeq2 (v1.42.1)
#   AggregateExpression -> DESeq2 default normalization
# ===========================================================================
rule pseudobulk_spatial:
    """
    Pseudo-bulk analysis for spatial data: aggregate raw counts per gene
    using Seurat AggregateExpression, then normalize with DESeq2 (v1.42.1)
    using default size factor estimation.
    """
    input:
        seurat_obj="results/spatial/seurat_spatial_clustered.rds",
    output:
        normalized="results/spatial/pseudobulk_deseq2_normalized.csv",
    log:
        "logs/pseudobulk_spatial.log",
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "scripts/pseudobulk_deseq2.R"


rule pseudobulk_snrna:
    """
    Pseudo-bulk analysis for snRNA-seq data: aggregate raw counts per gene
    using Seurat AggregateExpression, then normalize with DESeq2 (v1.42.1)
    using default size factor estimation.
    """
    input:
        seurat_obj="results/snrna/seurat_snrna_clustered.rds",
    output:
        normalized="results/snrna/pseudobulk_deseq2_normalized.csv",
    log:
        "logs/pseudobulk_snrna.log",
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "scripts/pseudobulk_deseq2.R"


# ===========================================================================
# STEP 11: hdWGCNA (v0.3.03) — co-expression network analysis
#   - Filter to soybean genes in >=1% of cells
#   - Subset to healthy + stressed mesophyll cells
#   - Signed network, minModuleSize=50, detectCutHeight=0.995,
#     mergeCutHeight=0.2
#   - TestSoftPowers for threshold selection
#   - Module eigengenes, hub genes, gene-centric UMAP
#   - Wilcoxon Rank Sum test for ME comparison
# ===========================================================================
rule hdwgcna_analysis:
    """
    High-dimensional weighted gene co-expression network analysis (hdWGCNA
    v0.3.03) on snRNA-seq data:
      1. Filter to soybean genes expressed in >=1% of cells
      2. Subset to healthy and stressed mesophyll cells
      3. TestSoftPowers to select soft-thresholding power
      4. ConstructNetwork: signed network, minModuleSize=50,
         detectCutHeight=0.995, mergeCutHeight=0.2
      5. Compute module eigengenes (MEs) and hub genes
      6. Gene-centric UMAP on topological overlap with top 10 hub genes
      7. Wilcoxon Rank Sum test for ME comparison between cell states
    """
    input:
        seurat_obj="results/snrna/seurat_snrna_clustered.rds",
    output:
        hdwgcna_obj="results/hdwgcna/hdwgcna_modules.rds",
        eigengenes="results/hdwgcna/module_eigengenes.csv",
        hub_genes="results/hdwgcna/hub_genes.csv",
    log:
        "logs/hdwgcna_analysis.log",
    threads: 4
    resources:
        mem_mb=32000,
    script:
        "scripts/hdwgcna_analysis.R"


# ===========================================================================
# STEP 12: Visualizations
#   - UMAP plots (Seurat)
#   - Heatmaps (pheatmap v1.0.12)
#   - Venn diagrams (eulerr v7.0.2 + SuperExactTest v1.1.0)
#   - UpSet plots (UpSetR v1.4.0)
#   - hdWGCNA gene-centric UMAP
# ===========================================================================
rule plot_umaps:
    """
    Generate UMAP visualizations for spatial and snRNA-seq data using Seurat
    DimPlot. Spatial UMAP colored by cluster and overlaid with ASR read
    percentage.
    """
    input:
        spatial_obj="results/spatial/seurat_spatial_clustered.rds",
        snrna_obj="results/snrna/seurat_snrna_clustered.rds",
    output:
        spatial_umap="results/figures/spatial_umap.pdf",
        snrna_umap="results/figures/snrna_umap.pdf",
    log:
        "logs/plot_umaps.log",
    script:
        "scripts/plot_umaps.R"


rule plot_marker_heatmap:
    """
    Heatmap of top cluster marker genes using pheatmap (v1.0.12).
    """
    input:
        spatial_markers="results/spatial/spatial_cluster_markers.csv",
        snrna_markers="results/snrna/snrna_cluster_markers.csv",
        spatial_obj="results/spatial/seurat_spatial_clustered.rds",
        snrna_obj="results/snrna/seurat_snrna_clustered.rds",
    output:
        heatmap="results/figures/marker_heatmap.pdf",
    log:
        "logs/plot_marker_heatmap.log",
    script:
        "scripts/plot_marker_heatmap.R"


rule plot_venn_and_upset:
    """
    Venn diagrams (eulerr v7.0.2 + SuperExactTest v1.1.0) and UpSet plots
    (UpSetR v1.4.0) comparing marker gene overlap between spatial clusters
    and snRNA-seq cell types.
    """
    input:
        spatial_markers="results/spatial/spatial_cluster_markers.csv",
        snrna_markers="results/snrna/snrna_cluster_markers.csv",
    output:
        venn="results/figures/venn_diagram.pdf",
        upset="results/figures/upset_plot.pdf",
    log:
        "logs/plot_venn_upset.log",
    script:
        "scripts/plot_venn_upset.R"


rule plot_hdwgcna:
    """
    Gene-centric UMAP of topological overlap with top 10 hub genes per
    module from hdWGCNA analysis.
    """
    input:
        hdwgcna_obj="results/hdwgcna/hdwgcna_modules.rds",
    output:
        gene_umap="results/figures/hdwgcna_gene_umap.pdf",
    log:
        "logs/plot_hdwgcna.log",
    script:
        "scripts/plot_hdwgcna.R"
