# Complete RNA-seq Analysis Tutorial for Local RStudio (Windows)

## System Requirements & Setup
- **RStudio Version**: 2025.05.1+513 "Mariposa Orchid"
- **System**: Windows 10/11, Intel Core i5
- **RAM**: Minimum 8GB recommended for this analysis
- **R Version**: 4.3.0 or higher recommended

## Dataset: Airway Smooth Muscle RNA-seq Data
We'll use the **airway** dataset - a gold standard RNA-seq dataset studying airway smooth muscle cells treated with dexamethasone. Perfect for learning because it's well-documented, moderate in size, and has clear biological interpretation.

---

## Step 1: Enhanced Package Installation for Local System

```r
# First, check your R version
R.version.string

# Set up CRAN mirror for faster downloads (choose your nearest mirror)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Check Bioconductor version
BiocManager::version()

# Core Bioconductor packages for RNA-seq analysis
bioc_packages <- c(
    "DESeq2",              # Main differential expression analysis
    "airway",              # Example dataset
    "org.Hs.eg.db",        # Human gene annotations
    "AnnotationDbi",       # Database interface
    "biomaRt",             # Gene annotations from Ensembl
    "GenomicFeatures",     # Genomic features manipulation
    "rtracklayer",         # Import/export genomic data
    "Rsamtools",           # BAM file handling
    "GenomicAlignments",   # Alignment data structures
    "SummarizedExperiment", # Data container
    "S4Vectors",           # S4 vector classes
    "IRanges",             # Integer ranges
    "GenomeInfoDb",        # Genome information
    "clusterProfiler",     # Gene ontology and pathway analysis
    "ReactomePA",          # Reactome pathway analysis
    "DOSE",                # Disease ontology analysis
    "enrichplot",          # Enhanced visualization for enrichment
    "pathview",            # Pathway visualization
    "ComplexHeatmap",      # Advanced heatmaps
    "EnhancedVolcano",     # Enhanced volcano plots
    "limma",               # Linear models for microarrays (useful utilities)
    "edgeR",               # Alternative DE analysis (for comparison)
    "tximport",            # Import transcript-level estimates
    "vsn"                  # Variance stabilizing normalization
)

# Install Bioconductor packages
BiocManager::install(bioc_packages, ask = FALSE, update = TRUE)

# Enhanced CRAN packages for better visualization and analysis
cran_packages <- c(
    "ggplot2",             # Grammar of graphics
    "dplyr",               # Data manipulation
    "tidyr",               # Data tidying
    "readr",               # Fast data reading
    "stringr",             # String manipulation
    "forcats",             # Factor manipulation
    "purrr",               # Functional programming
    "tibble",              # Modern data frames
    "plotly",              # Interactive plots
    "DT",                  # Interactive tables
    "corrplot",            # Correlation plots
    "VennDiagram",         # Venn diagrams
    "UpSetR",              # Upset plots for set intersections
    "pheatmap",            # Pretty heatmaps
    "RColorBrewer",        # Color palettes
    "viridis",             # Perceptually uniform color scales
    "scales",              # Scale functions for visualization
    "cowplot",             # Publication-ready plots
    "gridExtra",           # Grid arrangements
    "reshape2",            # Data reshaping
    "data.table",          # Fast data manipulation
    "openxlsx",            # Excel file handling
    "knitr",               # Dynamic reports
    "rmarkdown",           # R Markdown documents
    "flexdashboard",       # Interactive dashboards
    "shiny",               # Web applications
    "shinydashboard",      # Dashboard framework
    "ggsci",               # Scientific color palettes
    "ggrepel",             # Repel overlapping text labels
    "ggpubr",              # Publication-ready plots
    "dendextend",          # Dendrogram extensions
    "factoextra",          # Multivariate analysis visualization
    "cluster",             # Clustering algorithms
    "fpc",                 # Flexible procedures for clustering
    "NbClust",             # Optimal number of clusters
    "circlize",            # Circular visualization
    "wordcloud",           # Word clouds
    "treemap",             # Treemap visualization
    "networkD3",           # Network visualization
    "visNetwork",          # Network visualization
    "igraph",              # Graph analysis
    "survival",            # Survival analysis
    "survminer",           # Survival analysis visualization
    "broom",               # Tidy statistical output
    "modelr",              # Modeling functions
    "car",                 # Companion to applied regression
    "Hmisc",               # Miscellaneous functions
    "psych",               # Psychological research tools
    "GGally",              # Extension to ggplot2
    "corrr",               # Correlations in R
    "janitor",             # Data cleaning
    "skimr",               # Compact summary statistics
    "DataExplorer",        # Automated data exploration
    "mice",                # Multiple imputation
    "VIM",                 # Visualization of missing values
    "naniar",              # Explore missing data
    "here",                # Relative file paths
    "fs",                  # File system operations
    "glue",                # Glue strings
    "magrittr",            # Pipe operators
    "pryr",                # Pry into R
    "lobstr",              # Visualize R data structures
    "bench",               # Benchmarking
    "profvis",             # Profiling visualization
    "memoise",             # Memoization
    "future",              # Parallel processing
    "parallel",            # Parallel processing
    "doParallel",          # Parallel backend
    "foreach"              # Foreach looping
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Additional specialized packages for advanced analysis
advanced_packages <- c(
    "WGCNA",               # Weighted gene co-expression network analysis
    "GO.db",               # Gene Ontology database
    "KEGG.db",             # KEGG database
    "ReactomePA",          # Reactome pathway analysis
    "MSigDB",              # Molecular signatures database
    "fgsea",               # Fast gene set enrichment analysis
    "GSVA",                # Gene set variation analysis
    "sva",                 # Surrogate variable analysis
    "RUVSeq",              # Remove unwanted variation
    "NOISeq",              # Normalization and noise reduction
    "baySeq",              # Bayesian analysis of differential expression
    "EBSeq",               # Empirical Bayes analysis
    "ballgown",            # Analyze transcripts and genes
    "cummeRbund",          # Analysis of Cufflinks output
    "Glimma",              # Interactive plots
    "ReportingTools",      # Automated report generation
    "regioneR",            # Permutation tests for genomic regions
    "biovizBase",          # Visualization infrastructure
    "Gviz",                # Genomic data visualization
    "ggbio",               # Grammar of graphics for biological data
    "GenVisR",             # Genomic visualizations
    "maftools",            # Mutation annotation and visualization
    "TCGAbiolinks",        # TCGA data access
    "GEOquery",            # GEO data access
    "ArrayExpress",        # ArrayExpress data access
    "annotate",            # Annotation tools
    "GOstats",             # Gene ontology statistics
    "topGO",               # Topology-based GO analysis
    "Category",            # Category analysis
    "GSEABase",            # Gene set enrichment analysis base
    "graph",               # Graph data structures
    "Rgraphviz",           # Graph visualization
    "RBGL"                 # Boost graph library interface
)

# Install advanced packages (optional - install only if needed)
# BiocManager::install(advanced_packages, ask = FALSE)

# Check if all packages are installed successfully
required_packages <- c("DESeq2", "airway", "ggplot2", "dplyr", "pheatmap", "EnhancedVolcano")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) == 0) {
    cat("✓ All required packages are installed successfully!\n")
} else {
    cat("✗ Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Please install missing packages before proceeding.\n")
}

# Load core packages to verify installation
library(DESeq2)
library(airway)
library(ggplot2)
library(dplyr)

cat("✓ Core packages loaded successfully!\n")
cat("System ready for RNA-seq analysis!\n")
```

## Step 2: Set Up Working Environment

```r
# Create a dedicated project directory
project_dir <- "C:/Users/YourUsername/Documents/RNAseq_Analysis"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

# Create subdirectories for organized analysis
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("reports", showWarnings = FALSE)

# Set up parallel processing (utilize your Core i5 efficiently)
library(parallel)
library(BiocParallel)

# Detect number of cores
n_cores <- detectCores() - 1  # Leave one core for system
cat("Using", n_cores, "cores for parallel processing\n")

# Register parallel backend
register(MulticoreParam(workers = n_cores))

# Set up plotting parameters for high-quality output
options(repr.plot.width = 12, repr.plot.height = 8)
```

## Step 3: Load and Explore Data

```r
# Load required libraries
library(DESeq2)
library(airway)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(DT)
library(plotly)
library(ComplexHeatmap)
library(circlize)

# Load the airway dataset
data(airway)
se <- airway

# Enhanced data exploration
cat("=== DATASET OVERVIEW ===\n")
cat("Dataset dimensions:", dim(se), "\n")
cat("Number of samples:", ncol(se), "\n")
cat("Number of genes:", nrow(se), "\n")
cat("Memory usage:", format(object.size(se), units = "MB"), "\n\n")

# Sample information
cat("=== SAMPLE INFORMATION ===\n")
sample_info <- as.data.frame(colData(se))
print(sample_info)

# Treatment groups
cat("\n=== EXPERIMENTAL DESIGN ===\n")
cat("Treatment groups:\n")
print(table(se$dex))
cat("\nCell types:\n")
print(table(se$cell))
cat("\nDesign matrix:\n")
print(table(se$dex, se$cell))

# Gene information
cat("\n=== GENE INFORMATION ===\n")
cat("First few gene IDs:\n")
print(head(rownames(se), 10))

# Count matrix overview
cat("\n=== COUNT MATRIX PREVIEW ===\n")
print(assay(se)[1:6, 1:4])

# Basic statistics
cat("\n=== BASIC STATISTICS ===\n")
cat("Total counts per sample:\n")
print(colSums(assay(se)))
cat("\nMean counts per gene (top 10):\n")
print(head(sort(rowMeans(assay(se)), decreasing = TRUE), 10))
```

## Step 4: Enhanced Data Preprocessing

```r
# Create DESeq2 object with enhanced design
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Enhanced pre-filtering with multiple criteria
cat("=== PRE-FILTERING ANALYSIS ===\n")
cat("Genes before filtering:", nrow(dds), "\n")

# Filter 1: Remove genes with very low counts
keep_counts <- rowSums(counts(dds) >= 10) >= 3
cat("Genes passing count filter:", sum(keep_counts), "\n")

# Filter 2: Remove genes with low variance (optional)
rv <- rowVars(counts(dds))
keep_variance <- rv > quantile(rv, 0.1)
cat("Genes passing variance filter:", sum(keep_variance), "\n")

# Combine filters
keep_final <- keep_counts & keep_variance
dds <- dds[keep_final, ]

cat("Final genes after filtering:", nrow(dds), "\n")
cat("Percentage of genes retained:", round(nrow(dds)/nrow(se)*100, 1), "%\n")

# Set reference levels
dds$dex <- relevel(dds$dex, ref = "untrt")
dds$cell <- relevel(dds$cell, ref = "N61311")

# Save preprocessed data
save(dds, file = "data/preprocessed_dds.RData")
cat("✓ Preprocessed data saved\n")
```

## Step 5: Comprehensive Exploratory Data Analysis

```r
# Variance stabilizing transformation
cat("=== PERFORMING DATA TRANSFORMATION ===\n")
vst_data <- vst(dds, blind = FALSE)

# Enhanced PCA analysis
pca_data <- plotPCA(vst_data, intgroup = c("dex", "cell"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
p1 <- ggplot(pca_data, aes(PC1, PC2, color = dex, shape = cell)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = c("untrt" = "#E31A1C", "trt" = "#1F78B4")) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA Analysis: Treatment and Cell Type Effects",
         color = "Treatment", shape = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"))

print(p1)
ggsave("plots/enhanced_pca.png", p1, width = 10, height = 6, dpi = 300)

# Interactive PCA plot
p1_interactive <- ggplotly(p1, tooltip = c("colour", "shape"))
htmlwidgets::saveWidget(p1_interactive, "plots/interactive_pca.html")

# Sample correlation heatmap
sample_cor <- cor(assay(vst_data), method = "pearson")
pheatmap(sample_cor,
         annotation_col = as.data.frame(colData(vst_data)[, c("dex", "cell")]),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Sample Correlation Heatmap",
         filename = "plots/sample_correlation.png",
         width = 8, height = 6)

# Enhanced count distribution analysis
count_stats <- data.frame(
    Sample = colnames(dds),
    TotalCounts = colSums(counts(dds)),
    GenesDetected = colSums(counts(dds) > 0),
    MedianCounts = apply(counts(dds), 2, median),
    Treatment = dds$dex,
    CellType = dds$cell
)

# Count distribution plots
p2 <- ggplot(count_stats, aes(x = Sample, y = TotalCounts, fill = Treatment)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = c("untrt" = "#E31A1C", "trt" = "#1F78B4")) +
    labs(title = "Total Counts per Sample",
         x = "Sample", y = "Total Counts") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave("plots/count_distribution.png", p2, width = 10, height = 6, dpi = 300)

cat("✓ Exploratory analysis plots saved\n")
```

## Step 6: Differential Expression Analysis with Enhanced Statistics

```r
# Run DESeq2 with parallel processing
cat("=== RUNNING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
system.time({
    dds <- DESeq(dds, parallel = TRUE)
})

# Get results with multiple contrasts
# Primary contrast: treated vs untreated
res_primary <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
summary(res_primary)

# Enhanced results with shrinkage estimation
res_shrink <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

# Alternative shrinkage methods for comparison
res_normal <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "normal")
res_ashr <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "ashr")

# Compare shrinkage methods
par(mfrow = c(2, 2))
plotMA(res_primary, ylim = c(-5, 5), main = "Unshrunken")
plotMA(res_shrink, ylim = c(-5, 5), main = "apeglm")
plotMA(res_normal, ylim = c(-5, 5), main = "normal")
plotMA(res_ashr, ylim = c(-5, 5), main = "ashr")

# Use shrunken results for downstream analysis
res <- res_shrink
```

## Step 7: Enhanced Gene Annotation and Results Processing

```r
# Comprehensive gene annotation
cat("=== ADDING GENE ANNOTATIONS ===\n")

# Add multiple types of annotations
res$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res), 
                     column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$entrez <- mapIds(org.Hs.eg.db, keys = row.names(res), 
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
res$genename <- mapIds(org.Hs.eg.db, keys = row.names(res), 
                       column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
res$genetype <- mapIds(org.Hs.eg.db, keys = row.names(res), 
                       column = "GENETYPE", keytype = "ENSEMBL", multiVals = "first")

# Create comprehensive results table
results_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    dplyr::select(ensembl_id, symbol, genename, genetype, everything()) %>%
    arrange(padj) %>%
    filter(!is.na(padj))

# Enhanced filtering and categorization
results_df <- results_df %>%
    mutate(
        significance = case_when(
            padj < 0.001 & abs(log2FoldChange) > 2 ~ "Highly significant",
            padj < 0.01 & abs(log2FoldChange) > 1.5 ~ "Significant",
            padj < 0.05 & abs(log2FoldChange) > 1 ~ "Moderately significant",
            TRUE ~ "Not significant"
        ),
        regulation = case_when(
            log2FoldChange > 1 & padj < 0.05 ~ "Upregulated",
            log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
            TRUE ~ "Not regulated"
        ),
        effect_size = case_when(
            abs(log2FoldChange) > 2 ~ "Large",
            abs(log2FoldChange) > 1 ~ "Medium",
            abs(log2FoldChange) > 0.5 ~ "Small",
            TRUE ~ "Negligible"
        )
    )

# Summary statistics
cat("\n=== RESULTS SUMMARY ===\n")
cat("Total genes tested:", nrow(results_df), "\n")
cat("Significant genes (padj < 0.05):", sum(results_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated genes:", sum(results_df$regulation == "Upregulated"), "\n")
cat("Downregulated genes:", sum(results_df$regulation == "Downregulated"), "\n")

# Create interactive results table
datatable(results_df, 
          extensions = c('Buttons', 'ColReorder'),
          options = list(
              pageLength = 25,
              scrollX = TRUE,
              dom = 'Bfrtip',
              buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
              colReorder = TRUE,
              searchHighlight = TRUE
          ),
          caption = "Complete Differential Expression Results",
          filter = 'top') %>%
    formatRound(columns = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), 
                digits = 4) %>%
    formatStyle('padj', backgroundColor = styleInterval(0.05, c('lightcoral', 'white'))) %>%
    formatStyle('log2FoldChange', 
                backgroundColor = styleInterval(c(-1, 1), c('lightblue', 'white', 'lightcoral')))

# Save results
write.csv(results_df, "results/complete_results.csv", row.names = FALSE)
cat("✓ Results saved to complete_results.csv\n")
```

## Step 8: Advanced Visualization Suite

```r
# Enhanced volcano plot with multiple customizations
volcano_data <- results_df %>%
    filter(!is.na(padj)) %>%
    mutate(
        neg_log10_padj = -log10(padj),
        color_group = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
            padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
            TRUE ~ "Not significant"
        )
    )

# Create enhanced volcano plot
p_volcano <- EnhancedVolcano(volcano_data,
                            lab = volcano_data$symbol,
                            x = 'log2FoldChange',
                            y = 'padj',
                            selectLab = volcano_data$symbol[1:20],
                            xlab = bquote(~Log[2]~ 'fold change'),
                            ylab = bquote(~-Log[10]~adjusted~italic(P)),
                            title = 'Volcano Plot: Dexamethasone vs Untreated',
                            subtitle = 'Airway Smooth Muscle Cells',
                            pCutoff = 0.05,
                            FCcutoff = 1.0,
                            pointSize = 2.0,
                            labSize = 3.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colAlpha = 0.7,
                            legendPosition = 'right',
                            legendLabSize = 10,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            colConnectors = 'black',
                            max.overlaps = 20,
                            gridlines.major = FALSE,
                            gridlines.minor = FALSE,
                            border = 'partial',
                            borderWidth = 1.5,
                            borderColour = 'black')

print(p_volcano)
ggsave("plots/enhanced_volcano.png", p_volcano, width = 12, height = 8, dpi = 300)

# Advanced heatmap with ComplexHeatmap
top_genes <- results_df %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(50) %>%
    pull(ensembl_id)

heatmap_data <- assay(vst_data)[top_genes, ]
rownames(heatmap_data) <- results_df$symbol[match(rownames(heatmap_data), results_df$ensembl_id)]

# Create annotation
col_annotation <- HeatmapAnnotation(
    Treatment = dds$dex,
    CellType = dds$cell,
    col = list(
        Treatment = c("untrt" = "#E31A1C", "trt" = "#1F78B4"),
        CellType = c("N61311" = "#33A02C", "N052611" = "#FF7F00", 
                    "N080611" = "#6A3D9A", "N061011" = "#B15928")
    )
)

# Create complex heatmap
ht <- Heatmap(heatmap_data,
              name = "VST",
              top_annotation = col_annotation,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10),
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson",
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              heatmap_legend_param = list(title = "VST Expression"))

png("plots/complex_heatmap.png", width = 12, height = 10, units = "in", res = 300)
draw(ht)
dev.off()

cat("✓ Advanced visualizations saved\n")
```

## Step 9: Gene Set Enrichment Analysis

```r
# Prepare gene lists for enrichment analysis
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)

# Get significant genes
sig_genes_up <- results_df %>%
    filter(padj < 0.05, log2FoldChange > 1) %>%
    pull(entrez) %>%
    na.omit()

sig_genes_down <- results_df %>%
    filter(padj < 0.05, log2FoldChange < -1) %>%
    pull(entrez) %>%
    na.omit()

# Gene Ontology enrichment
ego_up <- enrichGO(gene = sig_genes_up,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = sig_genes_down,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)

# Visualize GO results
p_go_up <- barplot(ego_up, showCategory = 15, title = "GO Enrichment - Upregulated Genes")
p_go_down <- barplot(ego_down, showCategory = 15, title = "GO Enrichment - Downregulated Genes")

ggsave("plots/go_enrichment_up.png", p_go_up, width = 12, height = 8, dpi = 300)
ggsave("plots/go_enrichment_down.png", p_go_down, width = 12, height = 8, dpi = 300)

# KEGG pathway enrichment
kegg_up <- enrichKEGG(gene = sig_genes_up,
                      organism = "hsa",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")

kegg_down <- enrichKEGG(gene = sig_genes_down,
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")

# Save enrichment results
write.csv(as.data.frame(ego_up), "results/GO_enrichment_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down), "results/GO_enrichment_downregulated.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_up), "results/KEGG_enrichment_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_down), "results/KEGG_enrichment_downregulated.csv", row.names = FALSE)

cat("✓ Gene set enrichment analysis completed\n")
```

## Step 10: Generate Comprehensive Report

```r
# Create summary statistics
analysis_summary <- list(
    total_genes = nrow(results_df),
    significant_genes = sum(results_df$padj < 0.05, na.rm = TRUE),
    upregulated = sum(results_df$regulation == "Upregulated"),
    downregulated = sum(results_df$regulation == "Downregulated"),
    highly_significant = sum(results_df$significance == "Highly significant"),
    top_upregulated = results_df %>% filter(regulation == "Upregulated") %>% head(5) %>% pull(symbol),
    top_downregulated = results_df %>% filter(regulation == "Downregulated") %>% head(5) %>% pull(symbol)
)

# Print analysis summary
cat("=== FINAL ANALYSIS SUMMARY ===\n")
cat("Total genes analyzed:", analysis_summary$total_genes, "\n")
cat("Significant genes (padj < 0.05):", analysis_summary$significant_genes, "\n")
cat("Upregulated genes:", analysis_summary$upregulated, "\n")
cat("Downregulated genes:", analysis_summary$downregulated, "\n")
cat("Highly significant genes:", analysis_summary$highly_significant, "\n")
cat("\nTop 5 upregulated genes:", paste(analysis_summary$top_upregulated, collapse = ", "), "\n")
cat("Top 5 downregulated genes:", paste(analysis_summary$top_downregulated, collapse = ", "), "\n")

# Save session info and summary
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
save(analysis_summary, file = "results/analysis_summary.RData")
save.image("results/complete_analysis.RData")

# Export filtered results for different significance levels
highly_sig <- results_df %>% filter(padj < 0.001, abs(log2FoldChange) > 2)
moderate_sig <- results_df %>% filter(padj < 0.01, abs(log2FoldChange) > 1.5)
all_sig <- results_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)

write.csv(highly_sig, "results/highly_significant_genes.csv", row.names = FALSE)
write.csv(moderate_sig, "results/moderately_significant_genes.csv", row.names = FALSE)
write.csv(all_sig, "results/all_significant_genes.csv", row.names = FALSE)

cat("✓ Analysis complete! All results saved to respective folders.\n")
cat("✓ Check the 'results' folder for CSV files and 'plots' folder for visualizations.\n")
```

---

## Step 11: Performance Optimization & Memory Management

```r
# Monitor memory usage throughout analysis
cat("=== MEMORY USAGE MONITORING ===\n")
memory_usage <- function() {
    gc()
    cat("Memory usage:", format(object.size(ls(envir = .GlobalEnv)), units = "MB"), "\n")
    cat("R memory:", format(memory.size(), units = "MB"), "\n")
}

memory_usage()

# Optimize for your Intel Core i5 system
# Set up efficient parallel processing
library(parallel)
library(doParallel)

# Create cluster for parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Function to clean up large objects when not needed
cleanup_objects <- function() {
    # Remove large intermediate objects
    if(exists("vst_data")) rm(vst_data)
    if(exists("heatmap_data")) rm(heatmap_data)
    gc()
    cat("✓ Memory cleaned up\n")
}

# Stop cluster when done
# stopCluster(cl)
```

## Step 12: Quality Control and Validation

```r
# Additional quality control checks
cat("=== QUALITY CONTROL CHECKS ===\n")

# Check for batch effects
library(sva)
batch_effects <- ComBat_seq(counts(dds), batch = dds$cell, group = dds$dex)

# Validate results with alternative methods
library(edgeR)
library(limma)

# EdgeR analysis for comparison
y <- DGEList(counts = counts(dds), group = dds$dex)
y <- calcNormFactors(y)
design <- model.matrix(~ cell + dex, data = colData(dds))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "dextrt")
edger_results <- topTags(qlf, n = Inf)

# Compare DESeq2 and edgeR results
comparison_data <- merge(
    results_df[, c("ensembl_id", "symbol", "log2FoldChange", "padj")],
    edger_results$table[, c("logFC", "FDR")] %>% rownames_to_column("ensembl_id"),
    by = "ensembl_id", suffixes = c("_DESeq2", "_edgeR")
)

# Correlation between methods
cor_fc <- cor(comparison_data$log2FoldChange_DESeq2, comparison_data$logFC_edgeR, use = "complete.obs")
cor_p <- cor(-log10(comparison_data$padj_DESeq2), -log10(comparison_data$FDR_edgeR), use = "complete.obs")

cat("Fold change correlation between DESeq2 and edgeR:", round(cor_fc, 3), "\n")
cat("P-value correlation between DESeq2 and edgeR:", round(cor_p, 3), "\n")

# Create comparison plot
p_comparison <- ggplot(comparison_data, aes(x = log2FoldChange_DESeq2, y = logFC_edgeR)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = "DESeq2 vs edgeR: Fold Change Comparison",
         x = "DESeq2 log2FoldChange", y = "edgeR logFC") +
    theme_minimal()

ggsave("plots/method_comparison.png", p_comparison, width = 8, height = 6, dpi = 300)
```

## Step 13: Advanced Downstream Analysis

```r
# Time-course analysis simulation (if you have time-course data)
library(ImpulseDE2)
library(splines)

# Network analysis
library(WGCNA)
library(igraph)

# Perform WGCNA analysis
wgcna_data <- t(assay(vst_data))
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(wgcna_data, powerVector = powers, verbose = 5)

# Plot soft threshold
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")

# Construct network
net <- blockwiseModules(wgcna_data, power = 6, TOMType = "unsigned",
                       minModuleSize = 30, reassignThreshold = 0,
                       mergeCutHeight = 0.25, numericLabels = TRUE,
                       pamRespectsDendro = FALSE, saveTOMs = TRUE,
                       saveTOMFileBase = "results/TOM", verbose = 3)

# Plot module dendrogram
plotDendroAndColors(net$dendrograms[[1]], net$colors,
                   "Module colors", dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)

cat("✓ Network analysis completed\n")
```

## Step 14: Interactive Dashboard Creation

```r
# Create interactive dashboard
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)

# Dashboard UI
ui <- dashboardPage(
    dashboardHeader(title = "RNA-seq Analysis Dashboard"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
            menuItem("Results", tabName = "results", icon = icon("table")),
            menuItem("Visualization", tabName = "viz", icon = icon("chart-line")),
            menuItem("Enrichment", tabName = "enrichment", icon = icon("sitemap"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "overview",
                fluidRow(
                    box(title = "Analysis Summary", status = "primary", solidHeader = TRUE,
                        verbatimTextOutput("summary")),
                    box(title = "Sample Information", status = "info", solidHeader = TRUE,
                        DT::dataTableOutput("sample_info"))
                )
            ),
            tabItem(tabName = "results",
                fluidRow(
                    box(title = "Differential Expression Results", status = "primary", 
                        solidHeader = TRUE, width = 12,
                        DT::dataTableOutput("results_table"))
                )
            ),
            tabItem(tabName = "viz",
                fluidRow(
                    box(title = "Volcano Plot", status = "primary", solidHeader = TRUE,
                        plotlyOutput("volcano_plot")),
                    box(title = "MA Plot", status = "info", solidHeader = TRUE,
                        plotOutput("ma_plot"))
                )
            ),
            tabItem(tabName = "enrichment",
                fluidRow(
                    box(title = "GO Enrichment", status = "success", solidHeader = TRUE,
                        plotOutput("go_plot"))
                )
            )
        )
    )
)

# Save dashboard code
writeLines(capture.output(dput(ui)), "results/dashboard_ui.R")
cat("✓ Interactive dashboard code saved\n")
```

## Step 15: Automated Report Generation

```r
# Create automated R Markdown report
library(rmarkdown)
library(knitr)

# R Markdown template
rmd_template <- '
---
title: "RNA-seq Analysis Report: Airway Dataset"
author: "Bioinformatics Analysis"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
load("results/complete_analysis.RData")
```

## Executive Summary

This report presents the results of differential expression analysis on airway smooth muscle cells treated with dexamethasone.

### Key Findings:
- **Total genes analyzed**: `r analysis_summary$total_genes`
- **Significant genes**: `r analysis_summary$significant_genes`
- **Upregulated**: `r analysis_summary$upregulated`
- **Downregulated**: `r analysis_summary$downregulated`

## Sample Information

```{r sample-info}
knitr::kable(as.data.frame(colData(dds)), caption = "Sample Information")
```

## Results Overview

```{r results-summary}
sig_summary <- results_df %>%
  group_by(regulation) %>%
  summarise(Count = n()) %>%
  filter(regulation != "Not regulated")

ggplot(sig_summary, aes(x = regulation, y = Count, fill = regulation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "Significant Genes by Regulation Direction") +
  theme_minimal()
```

## Top Significant Genes

```{r top-genes}
top_genes_table <- results_df %>%
  filter(padj < 0.05) %>%
  select(symbol, genename, log2FoldChange, padj) %>%
  head(20)

knitr::kable(top_genes_table, digits = 4, caption = "Top 20 Significant Genes")
```

## Visualizations

```{r volcano-plot, fig.width=10, fig.height=6}
# Include volcano plot
knitr::include_graphics("plots/enhanced_volcano.png")
```

```{r heatmap, fig.width=10, fig.height=8}
# Include heatmap
knitr::include_graphics("plots/complex_heatmap.png")
```

## Gene Ontology Analysis

```{r go-analysis}
if(file.exists("results/GO_enrichment_upregulated.csv")) {
  go_up <- read.csv("results/GO_enrichment_upregulated.csv")
  knitr::kable(head(go_up, 10), caption = "Top GO Terms - Upregulated Genes")
}
```

## Session Information

```{r session-info}
sessionInfo()
```
'

# Write R Markdown file
writeLines(rmd_template, "reports/analysis_report.Rmd")

# Render the report
rmarkdown::render("reports/analysis_report.Rmd", 
                  output_file = "RNA_seq_Analysis_Report.html",
                  output_dir = "reports/")

cat("✓ Automated report generated: reports/RNA_seq_Analysis_Report.html\n")
```

## System-Specific Optimizations for Windows

```r
# Windows-specific optimizations
if(Sys.info()["sysname"] == "Windows") {
    # Set proper encoding
    options(encoding = "UTF-8")
    
    # Optimize for Windows file system
    options(stringsAsFactors = FALSE)
    
    # Set up proper paths for Windows
    if(!dir.exists("C:/RNAseq_temp")) {
        dir.create("C:/RNAseq_temp", showWarnings = FALSE)
    }
    
    # Set temporary directory
    Sys.setenv(TMPDIR = "C:/RNAseq_temp")
    
    # Memory optimization for Windows
    memory.limit(size = 8000)  # Set to 8GB if available
    
    cat("✓ Windows optimizations applied\n")
}
```

## Troubleshooting Guide

```r
# Common troubleshooting functions
troubleshoot <- function() {
    cat("=== TROUBLESHOOTING CHECKLIST ===\n")
    
    # Check R version
    cat("R version:", R.version.string, "\n")
    
    # Check memory
    cat("Available memory:", memory.limit(), "MB\n")
    
    # Check core packages
    required_packages <- c("DESeq2", "airway", "ggplot2", "dplyr")
    missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    
    if(length(missing) > 0) {
        cat("Missing packages:", paste(missing, collapse = ", "), "\n")
        cat("Install with: BiocManager::install(c('", paste(missing, collapse = "', '"), "'))\n")
    } else {
        cat("✓ All required packages available\n")
    }
    
    # Check file paths
    if(!dir.exists("results")) {
        dir.create("results", showWarnings = FALSE)
        cat("✓ Results directory created\n")
    }
    
    if(!dir.exists("plots")) {
        dir.create("plots", showWarnings = FALSE)
        cat("✓ Plots directory created\n")
    }
    
    cat("✓ Troubleshooting complete\n")
}

# Run troubleshooting
troubleshoot()
```

## Final Checklist and Next Steps

```r
# Final validation and next steps
cat("=== ANALYSIS COMPLETION CHECKLIST ===\n")

# Check all expected files exist
expected_files <- c(
    "results/complete_results.csv",
    "results/significant_genes.csv",
    "plots/enhanced_volcano.png",
    "plots/complex_heatmap.png",
    "results/analysis_summary.RData"
)

for(file in expected_files) {
    if(file.exists(file)) {
        cat("✓", file, "\n")
    } else {
        cat("✗", file, "(missing)\n")
    }
}

cat("\n=== NEXT STEPS FOR LEARNING ===\n")
cat("1. Explore the interactive results table in results/complete_results.csv\n")
cat("2. Examine the visualizations in the plots/ folder\n")
cat("3. Review the GO enrichment results\n")
cat("4. Try modifying parameters (p-value cutoffs, fold change thresholds)\n")
cat("5. Experiment with different gene sets for enrichment analysis\n")
cat("6. Practice with your own RNA-seq datasets\n")

cat("\n=== ADVANCED TOPICS TO EXPLORE ===\n")
cat("• Time-course RNA-seq analysis\n")
cat("• Single-cell RNA-seq analysis\n")
cat("• Multi-condition comparisons\n")
cat("• Pathway analysis and network reconstruction\n")
cat("• Integration with other omics data\n")

cat("\n✓ Tutorial completed successfully!\n")
cat("All results saved in:", getwd(), "\n")
```

---

## Performance Benchmarks for Your System

**Expected Runtime on Intel Core i5:**
- Package installation: 15-30 minutes (first time)
- Data preprocessing: 1-2 minutes
- DESeq2 analysis: 2-5 minutes
- Visualization generation: 3-5 minutes
- Complete analysis: 10-15 minutes

**Memory Usage:**
- Peak memory usage: ~2-4 GB
- Recommended minimum: 8 GB RAM
- Storage space needed: ~500 MB

**Optimization Tips:**
1. Use parallel processing (already included)
2. Save intermediate results frequently
3. Clean up large objects when not needed
4. Use efficient data structures
5. Consider using data.table for large datasets

This comprehensive tutorial is now optimized for your local Windows RStudio environment and will give you hands-on experience with professional-grade RNA-seq analysis workflows!