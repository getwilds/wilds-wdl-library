#!/usr/bin/env Rscript

library(scran)
library(SingleCellExperiment)

# Parse simple --key=value arguments (no optparse dependency; not present
# in the getwilds/scran image)
parse_args <- function(args, defaults) {
  opt <- defaults
  for (arg in args) {
    kv <- sub("^--", "", arg)
    key <- sub("=.*$", "", kv)
    value <- sub("^[^=]*=", "", kv)
    opt[[key]] <- value
  }
  opt
}

defaults <- list(
  sce_rds = NULL,
  sample_name = NULL,
  mito_pattern = "^MT-",
  nmads = "3.0",
  min_mean = "0.1",
  n_hvgs = "2000",
  output_prefix = NULL
)

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name
opt$nmads <- as.double(opt$nmads)
opt$min_mean <- as.double(opt$min_mean)
opt$n_hvgs <- as.integer(opt$n_hvgs)

set.seed(4)


##################
## 1. Load data ##
##################

message("Loading SingleCellExperiment from: ", opt$sce_rds)
sce <- readRDS(opt$sce_rds)
sce <- sce[!duplicated(rownames(sce)), ]
message("Cells loaded: ", ncol(sce))


##################
## 2. QC filter ##
##################

is_mito <- grepl(opt$mito_pattern, rownames(sce))
message("Mitochondrial genes detected: ", sum(is_mito))

counts_mat <- counts(sce)
lib_size <- colSums(counts_mat)
n_detected <- colSums(counts_mat > 0)
mito_percent <- if (sum(is_mito) > 0) {
  100 * colSums(counts_mat[is_mito, , drop = FALSE]) / lib_size
} else {
  rep(0, ncol(counts_mat))
}

# MAD-based outlier flagging on log-scale library size/detected genes, and
# raw-scale mitochondrial percentage
mad_outlier <- function(x, nmads, log = FALSE, side = "both") {
  vals <- if (log) log1p(x) else x
  med <- median(vals, na.rm = TRUE)
  spread <- mad(vals, center = med, na.rm = TRUE)
  lower <- med - nmads * spread
  upper <- med + nmads * spread
  switch(side,
    both = vals < lower | vals > upper,
    lower = vals < lower,
    higher = vals > upper
  )
}

discard <- mad_outlier(lib_size, opt$nmads, log = TRUE, side = "lower") |
  mad_outlier(n_detected, opt$nmads, log = TRUE, side = "lower") |
  mad_outlier(mito_percent, opt$nmads, side = "higher")
message("Cells flagged for discard: ", sum(discard), " / ", ncol(sce))

png(paste0(opt$output_prefix, "_qc.png"),
    width = 1400, height = 1200, res = 200)
plot(lib_size, n_detected,
     col = ifelse(discard, "firebrick", "grey40"),
     pch = 16, cex = 0.5,
     main = paste("Per-Cell QC:", opt$sample_name),
     xlab = "Total UMI Counts", ylab = "Detected Genes")
legend("bottomright", legend = c("kept", "discard"),
       col = c("grey40", "firebrick"), pch = 16, bty = "n")
dev.off()

sce <- sce[, !discard]
message("Cells after QC: ", ncol(sce))


#################################################
## 3. Normalize with scran deconvolution       ##
#################################################

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters, min.mean = opt$min_mean)
sce <- logNormCounts(sce)

png(paste0(opt$output_prefix, "_size_factors.png"),
    width = 1200, height = 1000, res = 200)
hist(sizeFactors(sce), breaks = 40,
     main = paste("Size Factor Distribution:", opt$sample_name),
     xlab = "Size Factor", ylab = "Number of Cells", col = "grey70")
dev.off()


#####################################################
## 4. Model variance and pick highly variable genes ##
#####################################################

dec <- modelGeneVar(sce, min.mean = opt$min_mean)
top_hvgs <- getTopHVGs(dec, n = opt$n_hvgs)

hvg_df <- as.data.frame(dec)
hvg_df$gene <- rownames(hvg_df)
write.csv(hvg_df[order(-hvg_df$bio), ],
          paste0(opt$output_prefix, "_hvgs.csv"),
          row.names = FALSE)

png(paste0(opt$output_prefix, "_mean_variance.png"),
    width = 1200, height = 1000, res = 200)
plot(hvg_df$mean, hvg_df$total,
     pch = 16, cex = 0.5, col = "grey40",
     main = paste("Mean-Variance Trend:", opt$sample_name),
     xlab = "Mean Log-Expression", ylab = "Variance")
points(hvg_df$mean[order(hvg_df$mean)], hvg_df$tech[order(hvg_df$mean)],
       type = "l", col = "dodgerblue", lwd = 2)
dev.off()


############################
## 5. PCA and save object ##
############################

# Base R PCA on the HVG subset of log-normalized counts (scran has no PCA
# function of its own; this keeps the script scran-only)
logcounts_hvg <- t(as.matrix(logcounts(sce)[top_hvgs, , drop = FALSE]))
pca_result <- prcomp(logcounts_hvg, center = TRUE, scale. = FALSE, rank. = 50)
reducedDims(sce) <- list(PCA = pca_result$x)

saveRDS(sce, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
