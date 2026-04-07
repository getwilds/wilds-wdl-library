version 1.0

import "ww-viennarna.wdl" as ww_viennarna
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow viennarna_example {
  # Generate test RNA sequences
  call ww_testdata.create_test_rna_fasta as download_demo_data { }

  # Predict MFE structures for all sequences in the FASTA
  call ww_viennarna.rnafold { input:
      input_fasta = download_demo_data.test_fasta,
      partition_function = true,
      cpu_cores = 1,
      memory_gb = 4
  }

  # Enumerate suboptimal structures for the test sequences
  call ww_viennarna.rnasubopt { input:
      input_fasta = download_demo_data.test_fasta,
      energy_range = 3.0,
      cpu_cores = 1,
      memory_gb = 4
  }

  # Predict RNA-RNA interaction between two short sequences
  call ww_viennarna.rnacofold { input:
      sequence_a = "GGGAAAUCCC",
      sequence_b = "GGGAUUUCCCC",
      cpu_cores = 1,
      memory_gb = 4
  }

  # Compute local base pair probabilities
  call ww_viennarna.rnaplfold { input:
      input_fasta = download_demo_data.test_fasta,
      cpu_cores = 1,
      memory_gb = 4
  }

  output {
    File rnafold_output = rnafold.structure_output
    Array[File] rnafold_plots = rnafold.postscript_plots
    File rnasubopt_output = rnasubopt.subopt_output
    File rnacofold_output = rnacofold.cofold_output
    Array[File] rnacofold_plots = rnacofold.postscript_plots
    Array[File] rnaplfold_profiles = rnaplfold.accessibility_profiles
    Array[File] rnaplfold_plots = rnaplfold.dp_plots
  }
}
