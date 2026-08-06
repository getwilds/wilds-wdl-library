#!/usr/bin/env Rscript

library(optparse)
library(scran)
library(scater)
library(DropletUtils)
library(SingleCellExperiment)
library(ggplot2)

option_list <- list(
  make_option("--input_h5",
              type = "character",
              default = NULL,
              help = "Path to Cell Ranger filtered feature-barcode matrix .h5 file"),
  make_option("--sample_name",
              type = "character",
              default = NULL,
              help = "Sample name for labeling outputs"),
  make_option("--mito_pattern",
              type = "character",
              default = "^MT-",
              help = "Regex pattern identifying mitochondrial gene symbols"),
  make_option("--nmads",
              type = "double",
              default = 3.0,
              help = "Number of MADs from the median used to flag low-quality cells"),
  make_option("--min_mean",
              type = "double",
              default = 0.1,
              help = "Minimum mean expression for a gene to be used in HVG/PCA selection"),
  make_option("--n_hvgs",
              type = "integer",
              default = 2000,
              help = "Number of top highly variable genes to select"),
  make_option("--output_prefix",
              type = "character",
              default = NULL,
              help = "Prefix for all output files (default: sample_name)")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name

set.seed(4)


##################
## 1. Load data ##
##################

message("Loading matrix from: ", opt$input_h5)
counts <- read10xCounts(opt$input_h5, col.names = TRUE)
counts <- counts[!duplicated(rownames(counts)), ]
message("Cells loaded: ", ncol(counts))


##################
## 2. QC filter ##
##################

is_mito <- grepl(opt$mito_pattern, rownames(counts))
message("Mitochondrial genes detected: ", sum(is_mito))

qc_stats <- perCellQCMetrics(counts, subsets = list(Mito = is_mito))
qc_filter <- quickPerCellQC(qc_stats,
                            sub.fields = "subsets_Mito_percent",
                            nmads = opt$nmads)

colData(counts) <- cbind(colData(counts), qc_stats)
counts$discard <- qc_filter$discard

# QC scatter plot of library size vs. detected features, colored by discard status
qc_plot <- plotColData(counts,
                       x = "sum",
                       y = "detected",
                       colour_by = "discard") +
  ggtitle(paste("Per-Cell QC for Sample:", opt$sample_name)) +
  xlab("Total UMI Counts") +
  ylab("Detected Genes")

ggsave(paste0(opt$output_prefix, "_qc.png"),
       plot = qc_plot,
       dpi = 300,
       width = 7,
       height = 6,
       device = "png")

sce <- counts[, !qc_filter$discard]
message("Cells after QC: ", ncol(sce))


#################################################
## 3. Normalize with scran deconvolution       ##
#################################################

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters, min.mean = opt$min_mean)
sce <- logNormCounts(sce)

sf_plot <- ggplot(data.frame(size_factor = sizeFactors(sce)),
                  aes(x = size_factor)) +
  geom_histogram(bins = 40) +
  ggtitle(paste("Size Factor Distribution for Sample:", opt$sample_name)) +
  xlab("Size Factor") +
  ylab("Number of Cells")

ggsave(paste0(opt$output_prefix, "_size_factors.png"),
       plot = sf_plot,
       dpi = 300,
       width = 6,
       height = 5,
       device = "png")


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

mean_var_plot <- ggplot(hvg_df, aes(x = mean, y = total)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = tech), color = "dodgerblue", linewidth = 1) +
  ggtitle(paste("Mean-Variance Trend for Sample:", opt$sample_name)) +
  xlab("Mean Log-Expression") +
  ylab("Variance")

ggsave(paste0(opt$output_prefix, "_mean_variance.png"),
       plot = mean_var_plot,
       dpi = 300,
       width = 6,
       height = 5,
       device = "png")


############################
## 5. PCA and save object ##
############################

sce <- runPCA(sce, subset_row = top_hvgs)

saveRDS(sce, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
