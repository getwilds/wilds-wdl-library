## WILDS WDL module for ViennaRNA
## Provides tasks for RNA secondary structure prediction and analysis using the
## ViennaRNA package, including minimum free energy folding, suboptimal structure
## enumeration, and RNA-RNA interaction prediction.

version 1.0

#### TASK DEFINITIONS ####

task rnafold {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Predict minimum free energy (MFE) secondary structures and partition function for RNA sequences using RNAfold"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-viennarna/ww-viennarna.wdl"
    outputs: {
        structure_output: "Text file containing predicted MFE structures in dot-bracket notation with free energies",
        postscript_plots: "Array of PostScript secondary structure plots (one per sequence)"
    }
    topic: "transcriptomics,nucleic_acid_structure_analysis"
    species: "human,eukaryote,prokaryote,virus"
    operation: "rna_secondary_structure_prediction"
    input_sample_required: "input_fasta:rna_sequence:fasta"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "structure_output:rna_secondary_structure:dot_bracket_format,postscript_plots:rna_secondary_structure_image:eps"
    output_reference: "none"
  }

  parameter_meta {
    input_fasta: "FASTA file containing RNA sequences to fold"
    partition_function: "Compute partition function and base pair probabilities in addition to MFE structure"
    no_lp: "Disallow lonely base pairs (isolated pairs with no adjacent pairs)"
    temperature: "Temperature in degrees Celsius for energy calculations"
    extra_args: "Additional command-line arguments to pass to RNAfold"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_fasta
    Boolean partition_function = false
    Boolean no_lp = false
    Float temperature = 37.0
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    RNAfold \
      --infile=~{input_fasta} \
      --temp=~{temperature} \
      ~{if partition_function then "--partfunc" else ""} \
      ~{if no_lp then "--noLP" else ""} \
      ~{extra_args} \
      > rnafold_output.txt
  >>>

  output {
    File structure_output = "rnafold_output.txt"
    Array[File] postscript_plots = glob("*_ss.ps")
  }

  runtime {
    docker: "getwilds/viennarna:2.7.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rnasubopt {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Enumerate suboptimal secondary structures within an energy range of the MFE using RNAsubopt"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-viennarna/ww-viennarna.wdl"
    outputs: {
        subopt_output: "Text file containing suboptimal structures in dot-bracket notation with energies"
    }
    topic: "transcriptomics,nucleic_acid_structure_analysis"
    species: "human,eukaryote,prokaryote,virus"
    operation: "rna_secondary_structure_prediction"
    input_sample_required: "input_fasta:rna_sequence:fasta"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "subopt_output:rna_secondary_structure:dot_bracket_format"
    output_reference: "none"
  }

  parameter_meta {
    input_fasta: "FASTA file containing an RNA sequence"
    energy_range: "Energy range above MFE in kcal/mol for suboptimal structure enumeration"
    no_lp: "Disallow lonely base pairs (isolated pairs with no adjacent pairs)"
    temperature: "Temperature in degrees Celsius for energy calculations"
    extra_args: "Additional command-line arguments to pass to RNAsubopt"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_fasta
    Float energy_range = 5.0
    Boolean no_lp = false
    Float temperature = 37.0
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    RNAsubopt \
      --infile=~{input_fasta} \
      --deltaEnergy=~{energy_range} \
      --temp=~{temperature} \
      ~{if no_lp then "--noLP" else ""} \
      ~{extra_args} \
      > rnasubopt_output.txt
  >>>

  output {
    File subopt_output = "rnasubopt_output.txt"
  }

  runtime {
    docker: "getwilds/viennarna:2.7.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rnacofold {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Predict the joint secondary structure and hybridization energy of two RNA molecules using RNAcofold"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-viennarna/ww-viennarna.wdl"
    outputs: {
        cofold_output: "Text file containing predicted joint MFE structure with interaction energy",
        postscript_plots: "Array of PostScript plots of the joint secondary structure"
    }
    topic: "transcriptomics,nucleic_acid_structure_analysis"
    species: "human,eukaryote,prokaryote,virus"
    operation: "rna_secondary_structure_prediction"
    input_sample_required: "none"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "cofold_output:rna_secondary_structure:dot_bracket_format,postscript_plots:rna_secondary_structure_image:eps"
    output_reference: "none"
  }

  parameter_meta {
    sequence_a: "First RNA sequence"
    sequence_b: "Second RNA sequence"
    partition_function: "Compute partition function in addition to MFE structure"
    temperature: "Temperature in degrees Celsius for energy calculations"
    extra_args: "Additional command-line arguments to pass to RNAcofold"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sequence_a
    String sequence_b
    Boolean partition_function = false
    Float temperature = 37.0
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # RNAcofold expects two sequences joined by '&'
    echo "~{sequence_a}&~{sequence_b}" | RNAcofold \
      --temp=~{temperature} \
      ~{if partition_function then "--partfunc" else ""} \
      ~{extra_args} \
      > rnacofold_output.txt
  >>>

  output {
    File cofold_output = "rnacofold_output.txt"
    Array[File] postscript_plots = glob("*_ss.ps")
  }

  runtime {
    docker: "getwilds/viennarna:2.7.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rnaplfold {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Compute local base pair probabilities for RNA sequences using a sliding window approach with RNAplfold"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-viennarna/ww-viennarna.wdl"
    outputs: {
        accessibility_profiles: "Array of accessibility profile files with unpaired probabilities",
        dp_plots: "Array of dot plot PostScript files with local pair probabilities"
    }
    topic: "transcriptomics,nucleic_acid_structure_analysis"
    species: "human,eukaryote,prokaryote,virus"
    operation: "rna_secondary_structure_prediction"
    input_sample_required: "input_fasta:rna_sequence:fasta"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "accessibility_profiles:rna_secondary_structure:textual_format,dp_plots:rna_secondary_structure_image:eps"
    output_reference: "none"
  }

  parameter_meta {
    input_fasta: "FASTA file containing RNA sequences"
    window_size: "Window size for the sliding window approach"
    span: "Maximum span (distance) of a base pair"
    unpaired_length: "Length of the unpaired region for accessibility computation"
    temperature: "Temperature in degrees Celsius for energy calculations"
    extra_args: "Additional command-line arguments to pass to RNAplfold"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_fasta
    Int window_size = 70
    Int span = 40
    Int unpaired_length = 4
    Float temperature = 37.0
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    RNAplfold \
      --winsize=~{window_size} \
      --span=~{span} \
      --ulength=~{unpaired_length} \
      --temp=~{temperature} \
      ~{extra_args} \
      < ~{input_fasta}
  >>>

  output {
    Array[File] accessibility_profiles = glob("*_lunp")
    Array[File] dp_plots = glob("*_dp.ps")
  }

  runtime {
    docker: "getwilds/viennarna:2.7.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
