version 1.0

# Import module under test using relative path for local development
import "ww-starling.wdl" as ww_starling

#### TEST WORKFLOW DEFINITION ####

workflow starling_example {
  # p53 N-terminal transactivation domain (residues 1-39) - a well-studied intrinsically disordered region
  String test_sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQA"

  # Generate ensemble for the test sequence
  call ww_starling.generate_ensemble { input:
    sequence = test_sequence,
    sample_name = "test_p53_ntad",
    num_conformations = 50,
    gpu_enabled = false,
    cpu_cores = 1,
    memory_gb = 4
  }

  # Query ensemble metadata
  call ww_starling.ensemble_info { input:
    starling_file = generate_ensemble.starling_file,
    sample_name = "test_p53_ntad"
  }

  output {
    File starling_file = generate_ensemble.starling_file
    File pdb_file = generate_ensemble.pdb_file
    File xtc_file = generate_ensemble.xtc_file
    File info_file = ensemble_info.info_file
  }
}
