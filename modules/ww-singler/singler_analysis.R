#!/usr/bin/env Rscript

library(scran)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# TEMPORARY: verbose HTTP tracing to debug an intermittent
# save_file()/array.h5 fetch failure seen only in Sprocket runs on HPC
# (not reproducible via GitHub CI or any manual repro attempt so far).
# Remove once root-caused.
library(httr)
httr::set_config(httr::verbose(data_out = TRUE, data_in = TRUE, info = TRUE, ssl = TRUE))
message("DEBUG HOME = ", Sys.getenv("HOME"))
writable <- function(p) if (nzchar(p)) file.access(p, mode = 2) == 0 else NA
for (v in c("GYPSUM_CACHE_DIR", "ANNOTATION_HUB_CACHE", "EXPERIMENT_HUB_CACHE")) {
  p <- Sys.getenv(v)
  message("DEBUG ", v, " = ", p,
          "  (exists=", dir.exists(p), " writable=", writable(p), ")")
}
message("DEBUG gypsum::cacheDirectory() = ",
        tryCatch(gypsum::cacheDirectory(), error = function(e) paste("ERR:", conditionMessage(e))))
message("DEBUG AnnotationHub cache dir = ",
        tryCatch(tools::R_user_dir("AnnotationHub", "cache"),
                 error = function(e) paste("ERR:", conditionMessage(e))))

# Parse simple --key=value arguments (no optparse dependency; not present
# in the getwilds/singler image)
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
  reference_rds = "",
  reference_dataset = "HumanPrimaryCellAtlasData",
  reference_ensembl = "true",
  label_column = "label.main",
  n_top_markers = "10",
  random_seed = "100",
  output_prefix = NULL
)

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name
opt$n_top_markers <- as.integer(opt$n_top_markers)
opt$random_seed <- as.integer(opt$random_seed)
opt$reference_ensembl <- as.logical(opt$reference_ensembl)

set.seed(opt$random_seed)


##################
## 1. Load data ##
##################

message("Loading SingleCellExperiment from: ", opt$sce_rds)
sce <- readRDS(opt$sce_rds)
message("Cells loaded: ", ncol(sce))

if (!"logcounts" %in% assayNames(sce)) {
  stop("Input SingleCellExperiment has no logcounts assay; ",
       "normalize (e.g. with ww-scran) before running singler")
}
if (!"PCA" %in% reducedDimNames(sce)) {
  stop("Input SingleCellExperiment has no PCA reducedDim; run ",
       "dimensionality reduction (e.g. with ww-scater) before running singler")
}


######################
## 2. Cluster cells ##
######################

clusters <- clusterCells(sce, use.dimred = "PCA")
sce$cluster <- clusters
message("Clusters identified: ", length(unique(clusters)))
write.csv(data.frame(barcode = colnames(sce), cluster = sce$cluster),
          paste0(opt$output_prefix, "_clusters.csv"),
          row.names = FALSE)


#################################
## 3. Per-cluster marker genes ##
#################################

markers <- findMarkers(sce, groups = sce$cluster)

# markers[[cluster]] is pre-sorted by ascending Top (best rank across all
# pairwise comparisons against other clusters); take the top N genes.
# Per-cluster tables have different logFC.<other cluster> columns, so only
# the shared summary columns are kept for the combined table
summary_cols <- c("Top", "p.value", "FDR", "summary.logFC")
marker_df <- do.call(rbind, lapply(names(markers), function(cluster_id) {
  top <- as.data.frame(markers[[cluster_id]])[, summary_cols]
  top <- head(top, opt$n_top_markers)
  data.frame(cluster = cluster_id, gene = rownames(top), top)
}))
write.csv(marker_df,
          paste0(opt$output_prefix, "_markers.csv"),
          row.names = FALSE)


##########################################
## 4. Cell-type annotation with SingleR ##
##########################################

if (nzchar(opt$reference_rds)) {
  message("Loading reference SingleCellExperiment from: ", opt$reference_rds)
  reference <- readRDS(opt$reference_rds)
} else {
  # 10x Cell Ranger output (the typical upstream source for sce) uses
  # Ensembl gene IDs as rownames, not gene symbols, so the celldex
  # reference is fetched in the matching ID space by default
  message("Fetching reference dataset from celldex: ", opt$reference_dataset,
          " (ensembl = ", opt$reference_ensembl, ")")
  reference <- do.call(opt$reference_dataset,
                        list(ensembl = opt$reference_ensembl))
}

predictions <- SingleR(test = sce, ref = reference,
                       labels = reference[[opt$label_column]],
                       clusters = sce$cluster)

sce$cell_type <- predictions$labels[match(sce$cluster, rownames(predictions))]

write.csv(as.data.frame(predictions),
          paste0(opt$output_prefix, "_predictions.csv"),
          row.names = TRUE)

saveRDS(sce, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
