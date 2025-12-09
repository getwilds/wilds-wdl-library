version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cnvkit/ww-cnvkit.wdl" as ww_cnvkit

workflow cnvkit_example {
  # Auto-download test data for testing purposes (multiple samples for reference)
  call ww_testdata.download_ref_data as download_ref { }
  call ww_testdata.download_bam_data as download_normal_1 {
    input: filename = "normal_sample_1.bam"
  }
  call ww_testdata.download_bam_data as download_normal_2 {
    input: filename = "normal_sample_2.bam"
  }
  call ww_testdata.download_bam_data as download_normal_3 {
    input: filename = "normal_sample_3.bam"
  }

  # Create reference using multiple normal samples
  call ww_cnvkit.create_reference {
    input:
      bam_files = [download_normal_1.bam, download_normal_2.bam, download_normal_3.bam],
      bam_indices = [download_normal_1.bai, download_normal_2.bai, download_normal_3.bai],
      reference_fasta = download_ref.fasta,
      reference_fasta_index = download_ref.fasta_index,
      cpu_cores = 1,
      memory_gb = 4
  }

  # Process single demo sample
  call ww_cnvkit.run_cnvkit {
    input:
      sample_name = "demo_sample",
      tumor_bam = download_normal_1.bam,
      tumor_bai = download_normal_1.bai,
      reference_cnn = create_reference.reference_cnn,
      cpu_cores = 1,
      memory_gb = 4
  }

  output {
    File cnv_calls = run_cnvkit.cnv_calls
    File cnv_segments = run_cnvkit.cnv_segments
    File cnv_plot = run_cnvkit.cnv_plot
    File reference_file = create_reference.reference_cnn
  }
}
