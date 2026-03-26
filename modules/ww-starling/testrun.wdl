version 1.0

# Import module under test and testdata module using relative paths for local development
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl" as ww_starling
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow starling_example {
  # p53 N-terminal transactivation domain (residues 1-39) - a well-studied intrinsically disordered region
  String test_sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQA"

  # Test generate_ensemble with a single sequence
  call ww_starling.generate_ensemble { input:
    sequence = test_sequence,
    sample_name = "test_p53_ntad",
    num_conformations = 50,
    gpu_enabled = false,
    cpu_cores = 1,
    memory_gb = 4
  }

  # Test ensemble_info on the generated ensemble
  call ww_starling.ensemble_info { input:
    starling_file = generate_ensemble.starling_file,
    sample_name = "test_p53_ntad"
  }

  # Test generate_ensemble_batch with a multi-sequence FASTA
  call ww_testdata.create_test_idp_fasta { }

  call ww_starling.generate_ensemble_batch { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    num_conformations = 50,
    gpu_enabled = false,
    cpu_cores = 1,
    memory_gb = 4
  }

  # Test split_fasta on the multi-sequence FASTA
  call ww_starling.split_fasta { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    sequences_per_batch = 2
  }

  output {
    File starling_file = generate_ensemble.starling_file
    File pdb_file = generate_ensemble.pdb_file
    File xtc_file = generate_ensemble.xtc_file
    File info_file = ensemble_info.info_file
    Array[File] batch_starling_files = generate_ensemble_batch.starling_files
    Array[File] batch_pdb_files = generate_ensemble_batch.pdb_files
    Array[File] batch_xtc_files = generate_ensemble_batch.xtc_files
    Array[File] split_batch_files = split_fasta.batch_files
  }
}
