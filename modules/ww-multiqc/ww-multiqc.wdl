## WILDS WDL MultiQC Module
## Aggregate bioinformatics analysis results into a single HTML report
## MultiQC searches directories for recognized log files from common bioinformatics
## tools and generates interactive summary visualizations across all samples.

version 1.0

#### TASK DEFINITIONS ####

task run_multiqc {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run MultiQC to aggregate QC reports from multiple bioinformatics tools into a single HTML report"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-multiqc/ww-multiqc.wdl"
    outputs: {
        html_report: "MultiQC interactive HTML summary report",
        data_dir: "Directory of parsed data files in tab-delimited format"
    }
    topic: "any,data_quality_management"
    species: "any"
    operation: "quality_control"
    in_sample_data: "quality_control_report"
    in_sample_format: "any"
    in_ref_data: "none"
    in_ref_format: "none"
    out_sample_data: "quality_control_report"
    out_sample_format: "html"
    out_ref_data: "none"
    out_ref_format: "none"
  }

  parameter_meta {
    input_files: "Array of files to scan for QC metrics (e.g., FastQC zips, alignment logs, etc.)"
    report_title: "Title displayed at the top of the MultiQC report"
    output_prefix: "Prefix for the output report filename"
    extra_args: "Additional command-line arguments to pass to MultiQC"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] input_files
    String report_title = "MultiQC Report"
    String output_prefix = "multiqc"
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Stage all input files into a single directory for MultiQC to scan
    mkdir -p qc_inputs
    for f in ~{sep=' ' input_files}; do
      ln -s "${f}" qc_inputs/
    done

    # Run MultiQC on the staged directory
    multiqc qc_inputs \
      --title "~{report_title}" \
      --filename "~{output_prefix}_report" \
      --outdir multiqc_output \
      --zip-data-dir \
      --force \
      ~{extra_args}

    echo "MultiQC report generated. Output files:"
    ls -la multiqc_output/
  >>>

  output {
    File html_report = "multiqc_output/~{output_prefix}_report.html"
    File data_dir = "multiqc_output/~{output_prefix}_report_data.zip"
  }

  runtime {
    docker: "getwilds/multiqc:1.33"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
