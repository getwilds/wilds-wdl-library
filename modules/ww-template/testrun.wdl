version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-template/ww-template.wdl" as ww_template
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structures for sample inputs if needed
struct TemplateSample {
    String name
    File input_file
}

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow template_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Create samples array using test data
  Array[TemplateSample] final_samples = [
    {
      "name": "demo_sample_1",
      "input_file": download_demo_data.r1_fastq
    },
    {
      "name": "demo_sample_2",
      "input_file": download_demo_data.r2_fastq
    }
  ]

  # Process each sample
  scatter (sample in final_samples) {
    call ww_template.process_sample { input:
        sample_name = sample.name,
        input_file = sample.input_file,
        cpu_cores = 1,
        memory_gb = 4
    }
  }

  output {
    Array[File] output_files = process_sample.output_file
  }
}
