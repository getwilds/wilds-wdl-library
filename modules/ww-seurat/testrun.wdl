version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-seurat/ww-seurat.wdl" as ww_seurat

workflow seurat_example {
  # TODO: Add test call of extract_h5

  # Download 10x Genomics H5 test data
  call ww_testdata.download_10x_h5_data { }

  call ww_seurat.run_seurat { input:
      input_h5 = download_10x_h5_data.h5_matrix,
      sample_name = "2500_Wistar_Rat_PBMCs_Singleplex"
  }

  output {
    File seurat_object = run_seurat.seurat_object
    File qc_plot = run_seurat.qc_plot
    File umap_plot = run_seurat.umap_plot
    File heatmap_plot = run_seurat.heatmap_plot
    File cluster_markers = run_seurat.cluster_markers
  }
}
