version 1.0

import "../../modules/ww-starling/ww-starling.wdl" as starling_tasks
import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow starling_batch_example {
  # Create a test FASTA with multiple short IDP sequences
  call ww_testdata.create_test_idp_fasta { }

  # Split into batches of 2 sequences each
  call starling_tasks.split_fasta { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    sequences_per_batch = 2
  }

  # Scatter ensemble generation across batches
  scatter (batch_fasta in split_fasta.batch_files) {
    call starling_tasks.generate_ensemble_batch { input:
      fasta_file = batch_fasta,
      num_conformations = 50,
      gpu_enabled = false,
      cpu_cores = 1,
      memory_gb = 4
    }
  }

  output {
    Array[Array[File]] starling_files = generate_ensemble_batch.starling_files
    Array[Array[File]] pdb_files = generate_ensemble_batch.pdb_files
    Array[Array[File]] xtc_files = generate_ensemble_batch.xtc_files
  }
}
