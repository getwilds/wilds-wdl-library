## WILDS WDL Template Module
## A simple functional template module that demonstrates WILDS WDL best practices.
## Contributors can copy this module and customize it for their specific tools.
## This template performs a simple "hello world" operation as a starting point.

version 1.0

# Import testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structures for sample inputs
struct TemplateSample {
    String name
    File input_file
}

#### WORKFLOW DEFINITION ####

workflow template_example {
  meta {
    author: "WILDS Development Team"
    email: "wilds@fredhutch.org"
    description: "WDL template module demonstrating WILDS best practices with simple hello world functionality"
    url: "https://github.com/getwilds/wilds-wdl-library"
    outputs: {
        output_files: "Array of simple output files with hello world message"
    }
  }

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
    call process_sample { input:
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

#### TASK DEFINITIONS ####

task process_sample {
  meta {
    description: "Simple template processing task that creates a hello world output file"
    outputs: {
        output_file: "Simple text file with hello world message and sample information"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_file: "Input file (any file type works for this template)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File input_file
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail
    
    # Simple hello world processing - replace this with your tool's commands
    echo "Hello, WILDS WDL World!" > "~{sample_name}.output.txt"
    echo "Sample: ~{sample_name}" >> "~{sample_name}.output.txt"
    echo "Input file: $(basename ~{input_file})" >> "~{sample_name}.output.txt"
    echo "Processing completed at: $(date)" >> "~{sample_name}.output.txt"
    
    # This is where you would put your actual tool commands, for example:
    # your_tool --input ~{input_file} --output ~{sample_name}.output.txt --threads ~{cpu_cores}
  >>>

  output {
    File output_file = "~{sample_name}.output.txt"
  }

  runtime {
    # Replace with your tool's Docker image
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

