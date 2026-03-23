version 1.0

import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "ww-starling-batch.wdl" as starling_batch_workflow

workflow starling_batch_example {
  # Create a test FASTA with multiple short IDP sequences
  call ww_testdata.create_test_idp_fasta { }

  # Run the batch pipeline
  call starling_batch_workflow.starling_batch { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    sequences_per_batch = 2,
    num_conformations = 50,
    gpu_enabled = false,
    cpu_cores = 1,
    memory_gb = 4
  }

  output {
    Array[Array[File]] starling_files = starling_batch.starling_files
    Array[Array[File]] pdb_files = starling_batch.pdb_files
    Array[Array[File]] xtc_files = starling_batch.xtc_files
  }
}
