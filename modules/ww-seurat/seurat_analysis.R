#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

option_list <- list(
  make_option("--input_h5",
              type = "character",
              default = NULL,
              help = "Path to Cell Ranger .h5 matrix file"),
  make_option("--sample_name",
              type = "character",
              default = NULL,
              help = "Sample name for labeling outputs"),
  make_option("--min_cells",
              type = "integer",
              default = 3,
              help = "Min cells per gene"),
  make_option("--min_features",
              type = "integer",
              default = 200,
              help = "Min features per cell"),
  make_option("--max_percent_mt",
              type = "double",
              default = 10.0,
              help = "Max mitochondrial percent"),
  make_option("--resolution",
              type = "double",
              default = 0.5,
              help = "Louvain clustering resolution"),
  make_option("--output_prefix",
              type = "character",
              default = NULL,
              help = "Prefix for all output files (default: sample_name)"),
  make_option("--ram_gb",
              type = "integer",
              default = 4,
              help = "Maximum RAM script can use in GB (default: 4)")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name
options(future.globals.maxSize = opt$ram_gb * 1024^3)

# Set max RAM usage (important for tasks like SCTransform and FindAllMarkers)
options(future.globals.maxSize = opt$ram_gb * 1000 * 1024^2)

set.seed(4)


#####################
## 1. Load H5 file ##
#####################

message("Loading matrix from: ", opt$input_h5)
counts <- Read10X_h5(opt$input_h5)

seurat_obj <- CreateSeuratObject(counts = counts,
                                 project = opt$sample_name,
                                 min.cells = opt$min_cells,
                                 min.features = opt$min_features)
message("Cells loaded: ", ncol(seurat_obj))


##################
## 2. QC filter ##
##################

# Get percent mitochondria (assume mitochondiral gene names start with "MT-")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Violin plot of QC metric
qc_plot <- VlnPlot(seurat_obj,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3,
                   pt.size  = 0.1) +
  plot_annotation(title = paste("Pre-filtering QC Metrics for Sample:", opt$sample_name))

# Save plot
ggsave(paste0(opt$output_prefix, "_qc.png"),
       plot = qc_plot,
       dpi = 300,
       width = 10,
       height = 5,
       device = "png")

# Remove cells outside of user-provided QC threholds
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > opt$min_features &
                              percent.mt < opt$max_percent_mt)
message("Cells after QC: ", ncol(seurat_obj))


#################################################
## 3. Normalize with SCTransform and make UMAP ##
#################################################

# Should have glmGamPoi installed for fast estimation
seurat_obj <- SCTransform(seurat_obj,
                          vars.to.regress = "percent.mt",
                          verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)


################
## 4. Cluster ##
################

seurat_obj <- FindClusters(seurat_obj,
                           resolution = opt$resolution,
                           verbose = FALSE)
message("Clusters found: ", length(unique(seurat_obj$seurat_clusters)))

# Make plot of UMAP colored by cluster
umap_plot <- DimPlot(seurat_obj,
                     pt.size = 0.4,
                     label = TRUE,
                     label.size = 5) +
  ggtitle(opt$sample_name) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2")

# Save plot
ggsave(paste0(opt$output_prefix, "_umap.png"),
       plot = umap_plot,
       dpi = 300,
       width = 6,
       height = 6,
       device = "png")


###########################################
## 5. Find marker genes and make heatmap ##
###########################################

message("Finding cluster marker genes...")
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          verbose = FALSE)

# Store top 30 markers per cluster to be saved for for manual annotation
top_30_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.csv(top_30_markers,
          paste0(opt$output_prefix, "_top30_markers.csv"),
          row.names = FALSE)

# Store top 8 markers for making the heatmap
top_8_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 8)

# Make a heatmap
heatmap_plot <- DoHeatmap(seurat_obj, features = top_8_markers$gene) +
  ggtitle(paste("Top 8 Markers Per Cluster for Sample:", opt$sample_name)) +
  NoLegend()

ggsave(paste0(opt$output_prefix, "_heatmap.png"),
       plot = heatmap_plot,
       dpi = 300,
       width = 10,
       height = 8,
       device = "png")


############################
## 6.  Save Seurat object ##
############################

saveRDS(seurat_obj, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
