version 1.0

import "ww-deeptools.wdl" as ww_deeptools
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

struct DeeptoolsSample {
    String name
    File bam
    File bai
}

#### TEST WORKFLOW DEFINITION ####

workflow deeptools_example {
  # Download reference data (includes BED file for regions)
  call ww_testdata.download_ref_data { }

  # Download two BAM files for multi-sample analysis
  call ww_testdata.download_bam_data as download_bam_1 { input:
      filename = "sample_1.bam"
  }
  call ww_testdata.download_bam_data as download_bam_2 { input:
      filename = "sample_2.bam"
  }

  # Create samples array
  Array[DeeptoolsSample] samples = [
    {
      "name": "sample_1",
      "bam": download_bam_1.bam,
      "bai": download_bam_1.bai
    },
    {
      "name": "sample_2",
      "bam": download_bam_2.bam,
      "bai": download_bam_2.bai
    }
  ]

  # Generate coverage tracks for each sample
  scatter (sample in samples) {
    call ww_deeptools.bam_coverage { input:
        bam = sample.bam,
        bai = sample.bai,
        sample_name = sample.name,
        normalize_using = "CPM",
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  # Compare two BAM files (using sample_1 as treatment, sample_2 as control)
  call ww_deeptools.bam_compare { input:
      treatment_bam = download_bam_1.bam,
      treatment_bai = download_bam_1.bai,
      control_bam = download_bam_2.bam,
      control_bai = download_bam_2.bai,
      sample_name = "test",
      normalize_using = "CPM",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Compute matrix of signal over genomic regions
  call ww_deeptools.compute_matrix { input:
      bigwig_files = bam_coverage.coverage_file,
      regions_file = download_ref_data.bed,
      sample_name = "test",
      mode = "scale-regions",
      before_region = 1000,
      after_region = 1000,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Generate heatmap from the matrix
  call ww_deeptools.plot_heatmap { input:
      matrix_gz = compute_matrix.matrix_gz,
      sample_name = "test",
      memory_gb = 4
  }

  # Generate profile plot from the matrix
  call ww_deeptools.plot_profile { input:
      matrix_gz = compute_matrix.matrix_gz,
      sample_name = "test",
      memory_gb = 4
  }

  # Multi-sample BAM summary for correlation and PCA
  call ww_deeptools.multi_bam_summary { input:
      bam_files = [download_bam_1.bam, download_bam_2.bam],
      bai_files = [download_bam_1.bai, download_bam_2.bai],
      sample_name = "test",
      bin_size = 50000,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Plot correlation from multi-sample summary
  call ww_deeptools.plot_correlation { input:
      summary_npz = multi_bam_summary.summary_npz,
      sample_name = "test",
      memory_gb = 4
  }

  # Plot PCA from multi-sample summary
  call ww_deeptools.plot_pca { input:
      summary_npz = multi_bam_summary.summary_npz,
      sample_name = "test",
      memory_gb = 4
  }

  # Plot fingerprint for QC
  call ww_deeptools.plot_fingerprint { input:
      bam_files = [download_bam_1.bam, download_bam_2.bam],
      bai_files = [download_bam_1.bai, download_bam_2.bai],
      sample_name = "test",
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    Array[File] coverage_files = bam_coverage.coverage_file
    File comparison = bam_compare.comparison_file
    File matrix = compute_matrix.matrix_gz
    File heatmap = plot_heatmap.heatmap
    File profile = plot_profile.profile
    File summary = multi_bam_summary.summary_npz
    File correlation = plot_correlation.correlation_plot
    File pca = plot_pca.pca_plot
    File fingerprint = plot_fingerprint.fingerprint_plot
  }
}
