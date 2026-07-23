#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
  library(RColorBrewer)
})

option_list <- list(
  make_option("--counts_file", type="character", help="Path to input counts matrix file (tab-delimited, genes as rows)"),
  make_option("--metadata_file", type="character", help="Path to sample metadata file (tab-delimited)"),
  make_option("--condition_column", type="character", default="condition", help="Column in metadata to use for comparison [default: %default]"),
  make_option("--reference_level", type="character", default="", help="Reference level for comparison [default: first alphabetically]"),
  make_option("--contrast", type="character", default="", help="Contrast: 'treatment,control' [default: infer from condition_column]"),
  make_option("--fdr_threshold", type="double", default=0.05, help="FDR threshold for significance [default: %default]"),
  make_option("--lfc_threshold", type="double", default=1.0, help="Minimum absolute log2 fold change for significance [default: %default]"),
  make_option("--min_count", type="integer", default=10, help="Minimum count for filterByExpr [default: %default]"),
  make_option("--min_total_count", type="integer", default=15, help="Minimum total count for filterByExpr [default: %default]"),
  make_option("--output_prefix", type="character", default="edger_results", help="Prefix for output files [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts_file) || is.null(opt$metadata_file)) {
  stop("Counts file and metadata file are required")
}

cat("Reading input data...\n")
counts_data <- read.delim(opt$counts_file, row.names=1, check.names=FALSE)
sample_metadata <- read.delim(opt$metadata_file, row.names=1, check.names=FALSE)

common_samples <- intersect(colnames(counts_data), rownames(sample_metadata))
if (length(common_samples) == 0) {
  stop("No matching samples between counts and metadata")
}
counts_data <- counts_data[, common_samples]
sample_metadata <- sample_metadata[common_samples, , drop=FALSE]

if (!opt$condition_column %in% colnames(sample_metadata)) {
  stop("Condition column '", opt$condition_column, "' not found in metadata")
}

condition <- as.factor(sample_metadata[[opt$condition_column]])
if (length(levels(condition)) < 2) {
  stop("Need at least two different levels in condition column for comparison")
}

if (opt$reference_level != "") {
  if (!opt$reference_level %in% levels(condition)) {
    stop("Reference level '", opt$reference_level, "' not found in condition column")
  }
  condition <- relevel(condition, ref=opt$reference_level)
}

cat("Setting up DGEList...\n")
dge <- DGEList(counts=round(counts_data), group=condition)

cat("Filtering low-count genes with filterByExpr...\n")
cat("Genes before filtering:", nrow(dge), "\n")
keep <- filterByExpr(dge, min.count=opt$min_count, min.total.count=opt$min_total_count)
dge <- dge[keep, , keep.lib.sizes=FALSE]
cat("Genes after filtering:", nrow(dge), "\n")

cat("Calculating TMM normalization factors...\n")
dge <- calcNormFactors(dge, method="TMM")

cat("Estimating dispersion...\n")
design <- model.matrix(~ condition)
dge <- estimateDisp(dge, design)

cat("Creating BCV plot...\n")
pdf(paste0(opt$output_prefix, "_bcv_plot.pdf"), width=8, height=6)
plotBCV(dge, main="Biological Coefficient of Variation")
dev.off()

cat("Fitting quasi-likelihood model...\n")
fit <- glmQLFit(dge, design)

cat("Running quasi-likelihood F-test...\n")
if (opt$contrast != "") {
  contrast_parts <- strsplit(opt$contrast, ",")[[1]]
  if (length(contrast_parts) != 2) {
    stop("Contrast should be in format: treatment,control")
  }
  treatment <- contrast_parts[1]
  control <- contrast_parts[2]
  contrast_levels <- levels(condition)
  treatment_idx <- which(contrast_levels == treatment)
  control_idx <- which(contrast_levels == control)
  if (length(treatment_idx) == 0 || length(control_idx) == 0) {
    stop("One or both contrast levels not found in condition column")
  }
  contrast_vec <- makeContrasts(
    contrasts=paste0("condition", treatment, "-condition", control),
    levels=design
  )
  qlf <- glmQLFTest(fit, contrast=contrast_vec)
} else {
  qlf <- glmQLFTest(fit, coef=2)
}

cat("Extracting results...\n")
all_results <- topTags(qlf, n=Inf, sort.by="PValue")$table
all_results$gene <- rownames(all_results)

sig_results <- subset(all_results, FDR < opt$fdr_threshold & abs(logFC) >= opt$lfc_threshold)
cat("Significant genes (FDR <", opt$fdr_threshold, ", |LFC| >=", opt$lfc_threshold, "):", nrow(sig_results), "\n")

cat("Computing normalized CPM...\n")
norm_cpm <- cpm(dge, normalized.lib.sizes=TRUE, log=FALSE)

cat("Writing output files...\n")
write.csv(all_results, file=paste0(opt$output_prefix, "_all_genes.csv"), row.names=FALSE)
write.csv(sig_results, file=paste0(opt$output_prefix, "_significant.csv"), row.names=FALSE)
write.csv(norm_cpm, file=paste0(opt$output_prefix, "_normalized_cpm.csv"))

cat("Creating MD plot...\n")
pdf(paste0(opt$output_prefix, "_md_plot.pdf"), width=8, height=6)
plotMD(qlf, main="MD Plot")
abline(h=c(-opt$lfc_threshold, opt$lfc_threshold), col="blue", lty=2)
dev.off()

cat("Creating volcano plot...\n")
volcano_data <- all_results
volcano_data$significant <- ifelse(
  volcano_data$FDR < opt$fdr_threshold & abs(volcano_data$logFC) >= opt$lfc_threshold,
  "Significant", "Not Sig"
)
volcano_data$neg_log10_fdr <- -log10(volcano_data$FDR)

volcano_plot <- ggplot(volcano_data, aes(x=logFC, y=neg_log10_fdr, color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("Significant"="red", "Not Sig"="grey")) +
  theme_classic() +
  geom_vline(xintercept=c(-opt$lfc_threshold, opt$lfc_threshold), linetype="dashed") +
  geom_hline(yintercept=-log10(opt$fdr_threshold), linetype="dashed") +
  labs(x="Log2 Fold Change", y="-Log10 FDR", color="Significance") +
  ggtitle("Volcano Plot")
ggsave(paste0(opt$output_prefix, "_volcano.pdf"), volcano_plot, width=8, height=6)

cat("Creating heatmap...\n")
if (nrow(sig_results) > 0) {
  top_genes <- rownames(sig_results)[1:min(50, nrow(sig_results))]
  log_cpm <- cpm(dge, normalized.lib.sizes=TRUE, log=TRUE, prior.count=2)
  heatmap_data <- log_cpm[top_genes, ]
  heatmap_data_z <- t(scale(t(heatmap_data)))

  anno_col <- data.frame(Condition=sample_metadata[[opt$condition_column]])
  rownames(anno_col) <- colnames(heatmap_data)

  pheatmap(
    heatmap_data_z,
    annotation_col=anno_col,
    clustering_distance_rows="correlation",
    clustering_distance_cols="correlation",
    show_rownames=TRUE,
    show_colnames=TRUE,
    fontsize_row=8,
    color=colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    filename=paste0(opt$output_prefix, "_heatmap.pdf"),
    width=8,
    height=10
  )
} else {
  cat("No significant genes found for heatmap.\n")
  pdf(paste0(opt$output_prefix, "_heatmap.pdf"), width=8, height=6)
  plot.new()
  text(0.5, 0.5, "No significant differentially expressed genes found")
  dev.off()
}

cat("edgeR analysis complete!\n")
