version 1.0

import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "./ww-pairtree-lineage.wdl" as pairtree_lineage_workflow

workflow pairtree_lineage_example {
  # Auto-generate two synthetic daughter-cell-line VCFs for testing purposes
  call ww_testdata.create_pairtree_vcfs { }

  # Run the clonal lineage inference pipeline
  call pairtree_lineage_workflow.pairtree_lineage { input:
    sample_vcfs = [create_pairtree_vcfs.sampleA_vcf, create_pairtree_vcfs.sampleB_vcf],
    sample_names = ["sampleA", "sampleB"],
    trees_per_chain = 100,
    cpu_cores = 1,
    memory_gb = 4
  }

  output {
    File ssm_file = pairtree_lineage.ssm_file
    File clustered_params_file = pairtree_lineage.clustered_params_file
    File results_file = pairtree_lineage.results_file
    File tree_html = pairtree_lineage.tree_html
    File tree_json = pairtree_lineage.tree_json
  }
}
