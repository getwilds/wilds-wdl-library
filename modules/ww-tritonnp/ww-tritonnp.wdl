## WILDS WDL for TritonNP - nucleosome positioning analysis from cfDNA
## Generates phasing features using FFT-based fragment size analysis

version 1.0

task triton_main {
  meta {
    author: "Chris Lo"
    email: "clo2@fredhutch.org"
    description: "Task for running TritonNP"
    outputs: {
        fm_file: "Array of outputs from TritonNP"
    }
  }

  parameter_meta {
    sample_name: "Sample name"
    bam_path: "BAM file"
    bam_index_path: "BAM index file"
    bias_path: "GC corrected file from Griffin"
    annotation: "BED file of genomic region to process on"
    reference_genome: "Reference genome file"
    reference_genome_index: "Reference genome file index"
    results_dir: "Output directory name"
    plot_list: "File containing names of genes to plot."
    map_quality: "Mapping quality threshold as a positive integer"
    size_range: "Size range as a space-delimited string, such as '15 500'"
    ncpus: "Number of CPUs to use"
    memory_gb: "Memory allocated for the task in GB"
  }
  
  input {
    String sample_name
    File bam_path
    File bam_index_path
    File bias_path
    File annotation
    File reference_genome
    File reference_genome_index
    String results_dir
    File plot_list
    Int map_quality = 20
    String size_range = "15 500"
    Int ncpus = 4
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Create results directory if it doesn't exist
    mkdir -p ~{results_dir}
    
    # Needed Python packages for running Triton
    python3 -m pip install pysam numpy scipy pandas matplotlib seaborn

    # Download script
    git clone https://github.com/caalo/TritonNP.git

    # Run Triton
    python TritonNP/GenerateFFTFeatures.py \
      --sample_name ~{sample_name} \
      --input ~{bam_path} \
      --bias ~{bias_path} \
      --annotation ~{annotation}  \
      --reference_genome ~{reference_genome} \
      --results_dir ~{results_dir} \
      --map_quality ~{map_quality} \
      --size_range ~{size_range} \
      --cpus ~{ncpus} \
      --plot_list ~{plot_list}
  >>>

  output {
    File fm_file = "~{results_dir}/~{sample_name}_PhasingFM.tsv"
  }

  runtime {
    cpu: ncpus
    memory: "~{memory_gb} GB"
    docker: "python:bullseye"
  }
}

task combine_fms {
  meta {
    author: "Chris Lo"
    email: "clo2@fredhutch.org"
    description: "Task for combine all sample outputs from TritonNP together"
    outputs: {
        final: "Aggregrated output file from TritonNP."
    }
  }

  parameter_meta {
    fm_files: "Array of output files from TritonNP"
    results_dir: "Output directory name"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] fm_files
    String results_dir
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Create results directory if it doesn't exist
    mkdir -p ~{results_dir}
    
    # Download script
    git clone https://github.com/caalo/TritonNP.git

    python3 -m pip install pandas

    # Run cleanup script
    python TritonNP/CombinePhasingFM.py \
       --inputs ~{sep=" " fm_files} \
       --results_dir ~{results_dir}
  >>>

  output {
    File final = "~{results_dir}/PhasingFM.tsv"
  }

  runtime {
    cpu: 1
    memory: "~{memory_gb} GB"
    docker: "python:bullseye"
  }
}
