version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "ww-clair3.wdl" as ww_clair3
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

struct Clair3Sample {
    String name
    File bam
    File bai
}

#### TEST WORKFLOW DEFINITION ####

workflow clair3_example {
  # Download reference genome and BAM test data
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1"
  }

  call ww_testdata.download_bam_data as download_bam_1 { input:
    filename = "sample1.bam"
  }

  call ww_testdata.download_bam_data as download_bam_2 { input:
    filename = "sample2.bam"
  }

  # Create samples array using test data
  Array[Clair3Sample] samples = [
    {
      "name": "sample1",
      "bam": download_bam_1.bam,
      "bai": download_bam_1.bai
    },
    {
      "name": "sample2",
      "bam": download_bam_2.bam,
      "bai": download_bam_2.bai
    }
  ]

  # Call Clair3 on each sample using Illumina model (test data is Illumina short-read)
  scatter (sample in samples) {
    call ww_clair3.run_clair3 { input:
      sample_name = sample.name,
      input_bam = sample.bam,
      input_bam_index = sample.bai,
      ref_fasta = download_reference.fasta,
      ref_fasta_index = download_reference.fasta_index,
      platform = "ilmn",
      model_path = "/opt/models/ilmn",
      ctg_name = "chr1",
      cpu_cores = 2,
      memory_gb = 8
    }
  }

  output {
    Array[File] output_vcfs = run_clair3.output_vcf
    Array[File] output_vcf_indices = run_clair3.output_vcf_index
  }
}
