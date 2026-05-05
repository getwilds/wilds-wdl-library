version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl" as sra_cellranger_workflow

#### TEST WORKFLOW DEFINITION ####
# HPC variant of the ww-sra-cellranger pipeline testrun. Dispatches via
# execution_mode = "hpc_sprocket" so the monthly HPC test run exercises
# the Sprocket + module-load path. Uses the same tiny SRA download and
# Cell Ranger reference as testrun.wdl so the only thing this run
# validates is the HPC dispatch path.

workflow sra_cellranger_example {
  call ww_testdata.download_test_cellranger_ref { }

  call sra_cellranger_workflow.sra_cellranger { input:
    sra_id_list = ["SRR9169219"],
    ref_gex = download_test_cellranger_ref.ref_tar,
    ncpu = 2,
    memory_gb = 6,
    max_reads = 100000,
    create_bam = false,
    chemistry = "SC3Pv2",
    execution_mode = "hpc_sprocket"
  }

  output {
    Array[File] cellranger_results = sra_cellranger.cellranger_results
    Array[File] cellranger_web_summaries = sra_cellranger.cellranger_web_summaries
    Array[File] cellranger_metrics = sra_cellranger.cellranger_metrics
  }
}
