version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-shapemapper/modules/ww-shapemapper/ww-shapemapper.wdl" as ww_shapemapper
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structures for sample inputs if needed
struct ShapeMapperSample {
    String name
    File target_fa
    File modified_r1
    File modified_r2
    File untreated_r1
    File untreated_r2
}

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow shapemapper_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Create samples array using test data
  # Note: For actual ShapeMapper analysis, you would need:
  # - A target RNA FASTA file
  # - Modified (chemically treated) sample FASTQ files (R1 and R2)
  # - Untreated control sample FASTQ files (R1 and R2)
  # This example uses demo FASTQ data as placeholders
  Array[ShapeMapperSample] final_samples = [
    {
      "name": "demo_sample_1",
      "target_fa": download_demo_data.r1_fastq,  # Placeholder - use actual target FASTA
      "modified_r1": download_demo_data.r1_fastq,
      "modified_r2": download_demo_data.r2_fastq,
      "untreated_r1": download_demo_data.r1_fastq,
      "untreated_r2": download_demo_data.r2_fastq
    }
  ]

  # Process each sample
  scatter (sample in final_samples) {
    call ww_shapemapper.run_shapemapper { input:
        sample_name = sample.name,
        target_fa = sample.target_fa,
        modified_r1 = sample.modified_r1,
        modified_r2 = sample.modified_r2,
        untreated_r1 = sample.untreated_r1,
        untreated_r2 = sample.untreated_r2,
        min_depth = 5000,
        is_amplicon = false
    }
  }

  output {
    Array[File] shape_files = run_shapemapper.shape_file
    Array[File] log_files = run_shapemapper.log_file
  }
}
