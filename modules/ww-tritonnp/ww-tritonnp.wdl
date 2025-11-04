

version 1.0

# Import testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/caalo/wilds-wdl-library/refs/heads/ww-TritonNP/modules/ww-testdata/ww-testdata.wdl" as ww_testdata


#### WORKFLOW DEFINITION ####

workflow tritonnp_example {
  meta {
    author: "Chris Lo"
    email: "clo2@fredhutch.org"
    description: "TritonNP module."
    url: "https://github.com/caalo/TritonNP"
    outputs: {
        final_composite: "Aggregrated outputs from all files",
        fm_files: "Array of outputs from individal samples"
    }
  }

  # Auto-download test data for testing purposes
  call ww_testdata.download_tritonnp_data as download_demo_data { }

  # Sample configuration
  Array[String] sample_names = ["NA12878"]
  Array[File] bam_paths = [download_demo_data.bam]
  Array[File] bam_index_paths = [download_demo_data.bam_index]
  Array[File] bias_paths = [download_demo_data.bias]
  
  # Global configuration parameters
  String results_dir = "results"
  File annotation = download_demo_data.annotation
  File reference_genome = download_demo_data.reference
  File reference_genome_index = download_demo_data.reference_index
  File plot_list = download_demo_data.plot_list
  Int map_quality = 20
  String size_range = "15 500"
  Int triton_main_ncpus = 4

  # Process each sample
  scatter (i in range(length(sample_names))) {
    call triton_main {
      input:
        sample_name = sample_names[i],
        bam_path = bam_paths[i],
        bam_index_path = bam_index_paths[i],
        bias_path = bias_paths[i],
        annotation = annotation,
        reference_genome = reference_genome,
        reference_genome_index = reference_genome_index,
        results_dir = results_dir,
        map_quality = map_quality,
        size_range = size_range,
        cpus = triton_main_ncpus,
        plot_list = plot_list
    }
  }

  # Combine all feature matrices after all samples are processed
  call combine_fms {
    input:
      fm_files = triton_main.fm_file,
      results_dir = results_dir
  }

  output {
    File final_composite = combine_fms.final
    Array[File] feature_files = triton_main.fm_file
  }
}

#### TASK DEFINITIONS ####

task triton_main {
  meta {
    description: "Task for running TritonNP"
    outputs: {
        fm_files: "Array of outputs from TritonNP"
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
    map_quality: "Mapping quality threshold as a positive integer"
    size_range: "Size range as a space-delimited string, such as '15 500'"
    cpus: "Number of CPUs to use"
    plot_list: "File containing names of genes to plot."
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
    Int map_quality
    String size_range
    Int cpus
    File plot_list
  }

  command <<<
    # Create results directory if it doesn't exist
    mkdir -p ~{results_dir}
    
    #Needed Python packges for running Triton
    python3 -m pip install pysam numpy scipy pandas matplotlib seaborn

    #download script
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
      --cpus ~{cpus} \
      --plot_list ~{plot_list}
  >>>

  output {
    File fm_file = "~{results_dir}/~{sample_name}_PhasingFM.tsv"
  }

  runtime {
    cpu: cpus
    memory: "4 GB" 
    docker: "python:bullseye" 
  }
}

task combine_fms {
  meta {
    description: "Task for combine all sample outputs from TritonNP together"
    outputs: {
        final: "Aggregrated output file from TritonNP."
    }
  }
  parameter_meta {
    fm_files: "Array of output files from TritonNP"
    results_dir: "Output directory name"
  }
  input {
    Array[File] fm_files
    String results_dir
  }

  command <<<
    # Create results directory if it doesn't exist
    mkdir -p ~{results_dir}
        
    #download script
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
    memory: "4 GB" 
    docker: "python:bullseye" 
  }
}

