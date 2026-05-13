version 1.0

# HPC variant of the ww-starling-batch pipeline testrun. Mirrors testrun.wdl
# but enables GPU on the underlying STARLING calls and uses the module's
# default 400 conformations — intended to validate end-to-end behavior on
# Fred Hutch HPC infrastructure (Sprocket or PROOF).

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-starling-batch/ww-starling-batch.wdl" as starling_batch_workflow

workflow starling_batch_example {
  call ww_testdata.create_test_idp_fasta { }

  call starling_batch_workflow.starling_batch { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    sequences_per_batch = 2,
    num_conformations = 400,
    gpu_enabled = true,
    cpu_cores = 4,
    memory_gb = 8
  }

  output {
    Array[Array[File]] starling_files = starling_batch.starling_files
    Array[Array[File]] pdb_files = starling_batch.pdb_files
    Array[Array[File]] xtc_files = starling_batch.xtc_files
  }
}
