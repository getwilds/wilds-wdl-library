## WILDS WDL module for UMI-tools.
## Provides tasks for deduplicating BAM files based on Unique Molecular Identifiers
## (UMIs) encoded in read names. Commonly used downstream of UMI-aware adapter
## trimmers (e.g. fastp --umi) for PRO-seq, RNA-seq, and other library types where
## PCR duplicates must be distinguished from true biological duplicates.

version 1.0

#### TASK DEFINITIONS ####

task dedup {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Deduplicates a coordinate-sorted, indexed BAM using umi_tools dedup. UMIs must already be encoded in the read name, separated from the read ID by the configured separator."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-umi-tools/ww-umi-tools.wdl"
    outputs: {
        deduped_bam: "Deduplicated BAM file",
        deduped_bai: "Index for the deduplicated BAM file",
        log: "umi_tools dedup log file (contains input/output read counts and per-position stats)"
    }
  }

  parameter_meta {
    input_bam: "Coordinate-sorted BAM file with UMIs encoded in read names"
    input_bai: "Index for the input BAM file"
    sample_name: "Sample name used for output file naming"
    paired: "If true, pass --paired to umi_tools dedup (required for paired-end data)"
    umi_separator: "Character separating the read ID from the UMI in the read name (default ':')"
    method: "umi_tools dedup method (unique, percentile, cluster, adjacency, directional)"
    extra_args: "Additional arguments to pass to umi_tools dedup"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_bam
    File input_bai
    String sample_name
    Boolean paired = true
    String umi_separator = ":"
    String method = "directional"
    String extra_args = ""
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # umi_tools requires the .bai to sit next to the .bam; stage both into the cwd.
    ln -s "~{input_bam}" "~{sample_name}.input.bam"
    ln -s "~{input_bai}" "~{sample_name}.input.bam.bai"

    umi_tools dedup \
      -I "~{sample_name}.input.bam" \
      -S "~{sample_name}.deduped.bam" \
      -L "~{sample_name}.dedup.log" \
      --umi-separator="~{umi_separator}" \
      --method="~{method}" \
      ~{if paired then "--paired" else ""} \
      ~{extra_args}

    samtools index -@ ~{cpu_cores} "~{sample_name}.deduped.bam"
  >>>

  output {
    File deduped_bam = "~{sample_name}.deduped.bam"
    File deduped_bai = "~{sample_name}.deduped.bam.bai"
    File log = "~{sample_name}.dedup.log"
  }

  runtime {
    docker: "getwilds/umitools:1.1.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
