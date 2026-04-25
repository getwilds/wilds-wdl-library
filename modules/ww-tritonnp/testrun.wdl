version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-tritonnp/ww-tritonnp.wdl" as ww_tritonnp

struct TritonSample {
    String name
    File bam
    File bam_index
    File bias
}

workflow tritonnp_example {
  meta {
    author: "Chris Lo"
    email: "clo2@fredhutch.org"
    description: "TritonNP module test workflow"
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

  # Create sample configuration array
  Array[TritonSample] samples = [
    object {
      name: "NA12878",
      bam: download_demo_data.bam,
      bam_index: download_demo_data.bam_index,
      bias: download_demo_data.bias
    }
  ]

  # Process each sample
  scatter (sample in samples) {
    call ww_tritonnp.triton_main {
      input:
        sample_name = sample.name,
        bam_path = sample.bam,
        bam_index_path = sample.bam_index,
        bias_path = sample.bias,
        annotation = download_demo_data.annotation,
        reference_genome = download_demo_data.reference,
        reference_genome_index = download_demo_data.reference_index,
        results_dir = "results",
        plot_list = download_demo_data.plot_list,
        ncpus = ncpus
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
