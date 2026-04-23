#!/usr/bin/env Rscript

# Script for running DESeq2 analysis from a WDL workflow
# This script handles:
# 1. Reading in count data and sample metadata
# 2. Setting up and running DESeq2 analysis
# 3. Outputting results files and visualizations

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
  library(RColorBrewer)
})

# Parse command line arguments
option_list <- list(
  make_option("--counts_file", type="character", help="Path to input counts matrix file"),
  make_option("--metadata_file", type="character", help="Path to sample metadata file"),
  make_option("--condition_column", type="character", default="condition", help="Column in metadata to use for comparison [default: %default]"),
  make_option("--reference_level", type="character", default="", help="Reference level for comparison [default: first alphabetically]"),
  make_option("--contrast", type="character", default="", help="Contrast to use (comma-separated: condition,treatment,control) [default: infer from condition_column]"),
  make_option("--output_prefix", type="character", default="deseq2_results", help="Prefix for output files [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$counts_file) || is.null(opt$metadata_file)) {
  stop("Counts file and metadata file are required")
}

# Read input data
cat("Reading input data...\n")
counts_data <- read.delim(opt$counts_file, row.names=1, check.names=FALSE)
sample_metadata <- read.delim(opt$metadata_file, row.names=1, check.names=FALSE)

# Ensure samples in counts match metadata
common_samples <- intersect(colnames(counts_data), rownames(sample_metadata))
if (length(common_samples) == 0) {
  stop("No matching samples between counts and metadata")
}

# Subset data to matching samples
counts_data <- counts_data[, common_samples]
sample_metadata <- sample_metadata[common_samples, , drop=FALSE]

# Ensure the condition column exists
if (!opt$condition_column %in% colnames(sample_metadata)) {
  stop("Condition column '", opt$condition_column, "' not found in metadata")
}

# Ensure condition column has at least two levels
if (length(unique(sample_metadata[[opt$condition_column]])) < 2) {
  stop("Need at least two different levels in condition column for comparison")
}

# Set the condition column as a factor
sample_metadata[[opt$condition_column]] <- as.factor(sample_metadata[[opt$condition_column]])

# Set reference level if specified
if (opt$reference_level != "") {
  if (!opt$reference_level %in% levels(sample_metadata[[opt$condition_column]])) {
    stop("Reference level '", opt$reference_level, "' not found in condition column")
  }
  sample_metadata[[opt$condition_column]] <- relevel(sample_metadata[[opt$condition_column]], ref=opt$reference_level)
}

# Create DESeq2 dataset
cat("Setting up DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_data), # Ensure counts are integers
  colData = sample_metadata,
  design = as.formula(paste("~", opt$condition_column))
)

# Filter out genes with too few counts
cat("Filtering low count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get results
cat("Extracting results...\n")
if (opt$contrast != "") {
  # Parse contrast
  contrast_parts <- strsplit(opt$contrast, ",")[[1]]
  if (length(contrast_parts) != 3) {
    stop("Contrast should be in format: condition,treatment,control")
  }
  contrast_use <- c(contrast_parts[1], contrast_parts[2], contrast_parts[3])
  res <- results(dds, contrast=contrast_use)
} else {
  # Use default contrast
  res <- results(dds)
}

# Add gene names to results
res$gene <- rownames(res)

# Sort by adjusted p-value
res_ordered <- res[order(res$padj),]

# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Write results to files
cat("Writing output files...\n")

# All genes
write.csv(as.data.frame(res_ordered), file=paste0(opt$output_prefix, "_all_genes.csv"), row.names=FALSE)

# Significantly differentially expressed genes
sig_genes <- subset(res_ordered, padj < 0.05)
write.csv(as.data.frame(sig_genes), file=paste0(opt$output_prefix, "_significant.csv"), row.names=FALSE)

# Normalized counts
write.csv(normalized_counts, file=paste0(opt$output_prefix, "_normalized_counts.csv"))

# Create PCA plot
cat("Creating PCA plot...\n")
# Use vst() for larger datasets (>1000 genes), rlog() for smaller datasets
# vst() is faster but requires sufficient genes; rlog() is more robust for small datasets
n_genes <- nrow(dds)
if (n_genes >= 1000) {
  cat("Using vst() transformation for", n_genes, "genes...\n")
  vsd <- vst(dds, blind=FALSE)
} else {
  cat("Using rlog() transformation for", n_genes, "genes (< 1000)...\n")
  vsd <- rlog(dds, blind=FALSE)
}
pcaData <- plotPCA(vsd, intgroup=opt$condition_column, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=get(opt$condition_column), shape=get(opt$condition_column))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color=opt$condition_column, shape=opt$condition_column) +
  theme_classic() +
  ggtitle("PCA Plot")
ggsave(paste0(opt$output_prefix, "_pca.pdf"), pca_plot, width=8, height=6)

# Create volcano plot
cat("Creating volcano plot...\n")
volcano_data <- as.data.frame(res)
volcano_data$significant <- ifelse(volcano_data$padj < 0.05, "FDR < 0.05", "Not Sig")
volcano_data$log10padj <- -log10(volcano_data$padj)

volcano_plot <- ggplot(volcano_data, aes(x=log2FoldChange, y=log10padj, color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("FDR < 0.05"="red", "Not Sig"="grey")) +
  theme_classic() +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted P-value", color="Significance") +
  ggtitle("Volcano Plot")
ggsave(paste0(opt$output_prefix, "_volcano.pdf"), volcano_plot, width=8, height=6)

# Create heatmap of top differentially expressed genes
if (nrow(sig_genes) > 0) {
  cat("Creating heatmap of top differentially expressed genes...\n")
  # Get top 50 genes or all significant genes if fewer
  top_genes <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
  heatmap_data <- assay(vsd)[top_genes, ]
  
  # Scale rows
  heatmap_data_z <- t(scale(t(heatmap_data)))
  
  # Create annotation for samples
  anno_col <- data.frame(Condition = sample_metadata[[opt$condition_column]])
  rownames(anno_col) <- colnames(heatmap_data)
  
  # Create heatmap
  pheatmap(
    heatmap_data_z,
    annotation_col = anno_col,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    filename = paste0(opt$output_prefix, "_heatmap.pdf"),
    width = 8,
    height = 10
  )
} else {
  cat("No significant genes found for heatmap.\n")
  # Create empty heatmap file
  pdf(paste0(opt$output_prefix, "_heatmap.pdf"), width=8, height=6)
  plot.new()
  text(0.5, 0.5, "No significant differentially expressed genes found")
  dev.off()
}

cat("DESeq2 analysis complete!\n")
