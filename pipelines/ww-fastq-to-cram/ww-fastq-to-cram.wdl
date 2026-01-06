version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gatk/ww-gatk.wdl" as gatk_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks

# Convert a group of paired fastq.gz files into an unmapped cram
# Uses the convention: READ_GROUP_NAME=~{sample_name}_~{group_name}

#### STRUCT DEFINITIONS

struct FastqGroup {
  String group_name # name identifier for this group of FASTQ files (e.g., flowcell, lane, or run name)
  Array[File] fastq_r1_locations # array of input R1 fastq file locations
  Array[File] fastq_r2_locations # array of input R2 fastq file locations
}

struct SampleData {
  String sample_name # sample name to insert into the read group header
  String? library_name # library name to place into the LB attribute in the read group header (optional, defaults to sample_name)
  String? sequencing_center # location where the sample was sequenced (optional, defaults to '.')
  Array[FastqGroup] fastq_groups # array of FASTQ file groups for this sample
}

#### WORKFLOW DEFINITION

workflow fastq_to_cram {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert paired FASTQ files to unmapped CRAM format using WILDS WDL modules"
    url: "https://github.com/getwilds/wilds-wdl-library/pipelines/ww-fastq-to-cram"
    outputs: {
        unmapped_crams: "Array of unmapped CRAM files for each sample",
        unmapped_cram_indexes: "Array of index files for each unmapped CRAM file",
        validation_reports: "Array of validation reports for each CRAM file"
    }
  }

  parameter_meta {
    batch_info: "Array of SampleData structs describing the relevant metadata for each sample"
    cpu_cores: "Number of CPU cores to use for processing"
    memory_gb: "Memory allocation in GB"
  }

  input {
    Array[SampleData] batch_info
    Int cpu_cores = 6
    Int memory_gb = 12
  }

  scatter (sample in batch_info) { # for every sample in your batch,
    String base_file_name = sample.sample_name

    scatter (group in sample.fastq_groups) { # and for every group of FASTQs for that sample,
      call gatk_tasks.fastq_to_bam { # take all the fastqs in that group and make an unmapped bam
        input:
          r1_fastq = group.fastq_r1_locations,
          r2_fastq = group.fastq_r2_locations,
          base_file_name = base_file_name + "_" + group.group_name,
          sample_name = sample.sample_name,
          library_name = sample.library_name,
          sequencing_center = sample.sequencing_center,
          read_group_name = sample.sample_name + "_" + group.group_name,
          cpu_cores = cpu_cores,
          memory_gb = memory_gb
      }
    } # End FASTQ group scatter

    call samtools_tasks.merge_bams_to_cram { # then for each sample, merge all unmapped bams into one unmapped cram
      input:
        bams_to_merge = fastq_to_bam.unmapped_bam,
        base_file_name = base_file_name,
        cpu_cores = cpu_cores,
        memory_gb = memory_gb
    }

    call gatk_tasks.validate_sam_file { # then validate to make sure it checks out
      input:
        input_file = merge_bams_to_cram.cram,
        base_file_name = base_file_name,
        mode = "SUMMARY",
        ignore_warnings = false
    }
  } # End sample scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] unmapped_crams = merge_bams_to_cram.cram
    Array[File] unmapped_cram_indexes = merge_bams_to_cram.crai
    Array[File] validation_reports = validate_sam_file.validation_report
  }
}
