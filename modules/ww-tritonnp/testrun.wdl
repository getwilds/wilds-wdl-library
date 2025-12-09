version 1.0

import "https://raw.githubusercontent.com/caalo/wilds-wdl-library/refs/heads/ww-TritonNP/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/caalo/wilds-wdl-library/refs/heads/ww-TritonNP/modules/ww-tritonnp/ww-tritonnp.wdl" as ww_tritonnp

workflow tritonnp_example {
  meta {
    author: "Chris Lo"
    email: "clo2@fredhutch.org"
    description: "TritonNP module."
    url: "https://github.com/caalo/TritonNP"
    outputs: {
        final_composite: "Aggregrated outputs from all files",
        fm_files: "Array of outputs from individal samples"
    }
  }

  input {
    Int ncpus = 4
  }

  # Auto-download test data for testing purposes
  call ww_testdata.download_tritonnp_data as download_demo_data { }

  # Sample configuration
  Array[String] sample_names = ["NA12878"]
  Array[File] bam_paths = [download_demo_data.bam]
  Array[File] bam_index_paths = [download_demo_data.bam_index]
  Array[File] bias_paths = [download_demo_data.bias]
  

  # Process each sample
  scatter (i in range(length(sample_names))) {
    call ww_tritonnp.triton_main {
      input:
        sample_name = sample_names[i],
        bam_path = bam_paths[i],
        bam_index_path = bam_index_paths[i],
        bias_path = bias_paths[i],
        annotation = download_demo_data.annotation,
        reference_genome = download_demo_data.reference,
        reference_genome_index = download_demo_data.reference_index,
        results_dir = "results",
        map_quality = 20,
        size_range = "15 500",
        cpus = ncpus,
        plot_list = download_demo_data.plot_list
    }
  }

  # Combine all feature matrices after all samples are processed
  call ww_tritonnp.combine_fms {
    input:
      fm_files = triton_main.fm_file,
      results_dir = "test_results"
  }

  output {
    File final_composite = combine_fms.final
    Array[File] fm_files = triton_main.fm_file
  }
}
