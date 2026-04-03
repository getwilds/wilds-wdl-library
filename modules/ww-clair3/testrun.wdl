version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-clair3/modules/ww-clair3/ww-clair3.wdl" as ww_clair3
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-clair3/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow clair3_example {
  # Download reference genome and BAM test data
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1"
  }

  call ww_testdata.download_bam_data as download_bam { input:
    filename = "sample1.bam"
  }

  # Call Clair3 using Illumina model (test data is Illumina short-read)
  call ww_clair3.run_clair3 { input:
    sample_name = "sample1",
    input_bam = download_bam.bam,
    input_bam_index = download_bam.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    platform = "ilmn",
    model_path = "/opt/models/ilmn",
    ctg_name = "chr1",
    pileup_only = true,
    cpu_cores = 2,
    memory_gb = 8
  }

  output {
    File output_vcf = run_clair3.output_vcf
    File output_vcf_index = run_clair3.output_vcf_index
  }
}
