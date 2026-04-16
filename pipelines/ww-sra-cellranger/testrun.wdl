version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-sra-cellranger/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl" as sra_cellranger_workflow

workflow sra_cellranger_example {
  # Download a small GEX reference for Cell Ranger
  call ww_testdata.download_test_cellranger_ref { }

  # Call the actual sra_cellranger workflow with test data
  # Using SRR7722937: Human Merkel cell carcinoma 10x Chromium 3' v2 scRNA-seq (GSE117988)
  # Standard read lengths (26bp R1, 98bp R2, 8bp index)
  # Limiting to 100k reads for fast testing while retaining enough for barcode detection
  call sra_cellranger_workflow.sra_cellranger { input:
    sra_id_list = ["SRR7722937"],
    ref_gex = download_test_cellranger_ref.ref_tar,
    ncpu = 2,
    memory_gb = 6,
    max_reads = 100000,
    create_bam = false,
    chemistry = "SC3Pv2"
  }

  output {
    Array[File] cellranger_results = sra_cellranger.cellranger_results
    Array[File] cellranger_web_summaries = sra_cellranger.cellranger_web_summaries
    Array[File] cellranger_metrics = sra_cellranger.cellranger_metrics
  }
}
