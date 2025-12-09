## WILDS WDL Template Module
## A simple functional template module that demonstrates WILDS WDL best practices.
## Contributors can copy this module and customize it for their specific tools.
## This template performs a simple "hello world" operation as a starting point.

# Specify WDL version, only 1.0 is supported currently
version 1.0

#### TASK DEFINITIONS ####
# Define tasks for each functionality of the tool represented by this module

task process_sample {
  # Provide metadata describing the purpose and authorship of the task
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Simple template processing task that creates a hello world output file"
    outputs: {
        output_file: "Simple text file with hello world message and sample information"
    }
  }

  # Provide parameter metadata for clarity on each input
  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_file: "Input file (any file type works for this template)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  # Specify inputs required for the task
  input {
    String sample_name
    File input_file
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  # Define the command section to execute the task's functionality
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

  # Define outputs produced by the task
  output {
    File output_file = "~{sample_name}.output.txt"
  }

  # Specify runtime requirements for the task
  runtime {
    # Replace with your tool's Docker image
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
