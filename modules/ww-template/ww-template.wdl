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
        output_files: "Array of simple output files with hello world message",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name and input files"
    cpu_cores: "Number of CPU cores allocated for each task"
    memory_gb: "Memory allocated for each task in GB"
  }

  input {
    Array[TemplateSample]? samples
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  # Auto-download test data if no samples provided
  if (!defined(samples)) {
    call ww_testdata.download_fastq_data as download_demo_data { }
  }

  # Create samples array - either from input or from test data download
  Array[TemplateSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample_1",
      "input_file": select_first([download_demo_data.r1_fastq])
    },
    {
      "name": "demo_sample_2",
      "input_file": select_first([download_demo_data.r2_fastq])
    }
  ]

  # Process each sample
  scatter (sample in final_samples) {
    call process_sample { input:
        sample_name = sample.name,
        input_file = sample.input_file,
        cpu_cores = cpu_cores,
        memory_gb = memory_gb
    }
  }

  # Validate all outputs
  call validate_outputs { input:
      output_files = process_sample.output_file
  }

  output {
    Array[File] output_files = process_sample.output_file
    File validation_report = validate_outputs.report
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

task validate_outputs {
  meta {
    description: "Simple validation to check that output files exist and are non-empty"
    outputs: {
        report: "Basic validation report confirming file existence"
    }
  }

  parameter_meta {
    output_files: "Array of output files to validate"
  }

  input {
    Array[File] output_files
  }

  command <<<
    set -eo pipefail
    
    {
      echo "WILDS Template Module Validation Report"
      echo "======================================="
      echo "Generated on: $(date)"
      echo ""
      echo "Samples processed: ~{length(output_files)}"
      echo ""
    } > validation_report.txt
    
    # Simple validation - check each file exists and is non-empty
    output_files=(~{sep=" " output_files})
    
    validation_passed=true
    
    for i in "${!output_files[@]}"; do
        output_file="${output_files[$i]}"
        
        echo "Checking output: $output_file" >> validation_report.txt
        
        if [[ -f "$output_file" && -s "$output_file" ]]; then
            file_size=$(wc -c < "$output_file")
            echo "  Output file exists and is non-empty ($file_size bytes)" >> validation_report.txt
        else
            echo "  Output file is missing or empty" >> validation_report.txt
            validation_passed=false
        fi
    done
    
    echo "" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
        echo "ALL VALIDATIONS PASSED" >> validation_report.txt
    else
        echo "SOME VALIDATIONS FAILED" >> validation_report.txt
        exit 1
    fi
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    # Lightweight Docker image for simple bash operations
    # Replace with a more appropriate image if your validation requires specific tools
    docker: "getwilds/bwa:0.7.17"
    cpu: 1
    memory: "2 GB"
  }
}
