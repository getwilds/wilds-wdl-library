version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-fastp/ww-fastp.wdl" as fastp_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-bowtie2/ww-bowtie2.wdl" as bowtie2_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-umi-tools/ww-umi-tools.wdl" as umi_tools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-deeptools/ww-deeptools.wdl" as deeptools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-multiqc/ww-multiqc.wdl" as multiqc_tasks

struct ProseqSample {
    String name
    File r1
    File r2
}

struct ProseqReferences {
    File experimental_fasta
    File spikein_merged_fasta
    File rdna_fasta
    String spikein_chrom_prefix
}

workflow proseq {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Paired-end PRO-seq pipeline with UMI deduplication and spike-in normalization, adapted from JAJ256/PROseq_alignment.sh. Covers UMI-aware trimming, rRNA depletion + spike-in + experimental alignment, UMI dedup, strand-specific 3' bigWigs, and aggregated QC."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-proseq/ww-proseq.wdl"
    outputs: {
        fastp_r1_trimmed: "UMI-extracted, adapter-trimmed R1 FASTQ files for each sample",
        fastp_r2_trimmed: "UMI-extracted, adapter-trimmed R2 FASTQ files for each sample",
        fastp_html: "fastp HTML reports for each sample",
        fastp_json: "fastp JSON reports for each sample",
        rrna_unaligned_r1: "R1 reads that failed to align to the rDNA reference (rRNA-depleted)",
        rrna_unaligned_r2: "R2 reads that failed to align to the rDNA reference (rRNA-depleted)",
        experimental_bam_prededup: "Pre-dedup BAM aligned to the experimental genome",
        experimental_bam_dedup: "UMI-deduplicated BAM aligned to the experimental genome",
        experimental_bai_dedup: "Index for the UMI-deduplicated experimental BAM",
        experimental_dedup_log: "umi_tools dedup log for the experimental BAM (input/output read counts and PCR-duplicate stats)",
        spikein_bam_prededup: "Pre-dedup BAM filtered to spike-in chromosomes",
        spikein_bam_dedup: "UMI-deduplicated spike-in BAM (used to derive normalization factors across samples)",
        spikein_bai_dedup: "Index for the UMI-deduplicated spike-in BAM",
        spikein_dedup_log: "umi_tools dedup log for the spike-in BAM",
        bigwig_forward: "Strand-specific (forward) bigWig of 3' read ends at single-base resolution",
        bigwig_reverse: "Strand-specific (reverse) bigWig of 3' read ends at single-base resolution",
        multiqc_report: "MultiQC interactive HTML report aggregating fastp metrics",
        multiqc_data: "MultiQC data directory"
    }
  }

  parameter_meta {
    samples: "List of paired-end PRO-seq samples (name + R1 + R2)"
    references: "PRO-seq reference set: experimental genome FASTA, spike-in-merged FASTA, rDNA FASTA, and the chromosome-name prefix used to flag spike-in contigs in the merged reference"
    umi_loc: "fastp UMI location. 'per_read', 'read1', 'read2', or '' to disable UMI extraction (default: 'per_read')."
    umi_len: "UMI length in basepairs (default: 6)"
    mapq_threshold: "Minimum MAPQ score for retained alignments"
    align_cpu: "CPU cores allocated to bowtie2 alignment tasks"
    align_memory_gb: "Memory (GB) allocated to bowtie2 alignment tasks"
  }

  input {
    Array[ProseqSample] samples
    ProseqReferences references
    String umi_loc = "per_read"
    Int umi_len = 6
    Int mapq_threshold = 10
    Int align_cpu = 4
    Int align_memory_gb = 8
  }

  # Build the three bowtie2 indices once, in parallel, before the per-sample scatter.
  call bowtie2_tasks.bowtie2_build as build_experimental { input:
      reference_fasta = references.experimental_fasta,
      index_prefix = "experimental",
      cpu_cores = align_cpu,
      memory_gb = align_memory_gb
  }
  call bowtie2_tasks.bowtie2_build as build_spikein { input:
      reference_fasta = references.spikein_merged_fasta,
      index_prefix = "spikein",
      cpu_cores = align_cpu,
      memory_gb = align_memory_gb
  }
  call bowtie2_tasks.bowtie2_build as build_rdna { input:
      reference_fasta = references.rdna_fasta,
      index_prefix = "rdna",
      cpu_cores = align_cpu,
      memory_gb = align_memory_gb
  }

  scatter (sample in samples) {
    # Step 1 — UMI-aware adapter trimming. fastp moves the UMI from the read sequence
    # into the read name (separated by ':'), the format umi_tools dedup consumes later.
    call fastp_tasks.fastp_paired { input:
        sample_name = sample.name,
        r1_fastq = sample.r1,
        r2_fastq = sample.r2,
        umi_loc = umi_loc,
        umi_len = umi_len
    }

    # Step 2 — rRNA depletion. Align to the rDNA reference; reads that *fail* to align
    # are passed to subsequent alignment steps. The mapped BAM is discarded.
    call bowtie2_tasks.bowtie2_align as deplete_rrna { input:
        bowtie2_index_tar = build_rdna.bowtie2_index_tar,
        index_prefix = "rdna",
        reads = fastp_paired.r1_trimmed,
        mates = fastp_paired.r2_trimmed,
        name = sample.name + ".rrna",
        preset = "fast-local",
        capture_unaligned = true,
        cpu_cores = align_cpu,
        memory_gb = align_memory_gb
    }

    # Step 3a — Spike-in alignment. Aligns rRNA-depleted reads to a merged
    # experimental+spike-in reference, then keeps only proper pairs above the MAPQ
    # threshold. Concordance flags match the original PROseq_alignment.sh script.
    call bowtie2_tasks.bowtie2_align as align_spikein { input:
        bowtie2_index_tar = build_spikein.bowtie2_index_tar,
        index_prefix = "spikein",
        reads = select_first([deplete_rrna.unaligned_r1]),
        mates = select_first([deplete_rrna.unaligned_r2]),
        name = sample.name + ".spikein_unfiltered",
        preset = "very-sensitive-local",
        min_mapq = mapq_threshold,
        samtools_filter_flags = "-f 2",
        extra_bowtie2_args = "--no-mixed --no-discordant",
        cpu_cores = align_cpu,
        memory_gb = align_memory_gb
    }

    # Step 3b — Filter the merged-reference BAM down to the spike-in contigs only.
    call samtools_tasks.filter_bam_by_chrom_prefix as keep_spikein_only { input:
        input_bam = align_spikein.sorted_bam,
        input_bai = align_spikein.sorted_bai,
        chrom_prefix = references.spikein_chrom_prefix,
        sample_name = sample.name + ".spikein"
    }

    # Step 4 — Experimental genome alignment with proper-pair + MAPQ filter.
    call bowtie2_tasks.bowtie2_align as align_experimental { input:
        bowtie2_index_tar = build_experimental.bowtie2_index_tar,
        index_prefix = "experimental",
        reads = select_first([deplete_rrna.unaligned_r1]),
        mates = select_first([deplete_rrna.unaligned_r2]),
        name = sample.name + ".experimental",
        preset = "sensitive-local",
        min_mapq = mapq_threshold,
        samtools_filter_flags = "-f 2",
        cpu_cores = align_cpu,
        memory_gb = align_memory_gb
    }

    # Step 5 — UMI-aware PCR deduplication on both BAMs.
    call umi_tools_tasks.dedup as dedup_experimental { input:
        input_bam = align_experimental.sorted_bam,
        input_bai = align_experimental.sorted_bai,
        sample_name = sample.name + ".experimental"
    }
    call umi_tools_tasks.dedup as dedup_spikein { input:
        input_bam = keep_spikein_only.filtered_bam,
        input_bai = keep_spikein_only.filtered_bai,
        sample_name = sample.name + ".spikein"
    }

    # Step 6 — Strand-specific 3' single-base bigWigs from the dedup'd experimental BAM.
    # --Offset 1 counts the 3' end of the read (the active site of Pol II in PRO-seq).
    # --samFlagInclude 82 selects forward-strand fragments, 98 selects reverse-strand.
    # --normalizeUsing None matches the script; downstream analyses can re-normalize
    # using spike-in-derived factors.
    call deeptools_tasks.bam_coverage as bigwig_fwd { input:
        bam = dedup_experimental.deduped_bam,
        bai = dedup_experimental.deduped_bai,
        sample_name = sample.name + ".fwd",
        bin_size = 1,
        normalize_using = "None",
        extra_args = "--Offset 1 --samFlagInclude 82 --skipNonCoveredRegions"
    }
    call deeptools_tasks.bam_coverage as bigwig_rev { input:
        bam = dedup_experimental.deduped_bam,
        bai = dedup_experimental.deduped_bai,
        sample_name = sample.name + ".rev",
        bin_size = 1,
        normalize_using = "None",
        extra_args = "--Offset 1 --samFlagInclude 98 --skipNonCoveredRegions"
    }
  }

  # Aggregate fastp metrics into one MultiQC report.
  call multiqc_tasks.run_multiqc { input:
      input_files = flatten([
        fastp_paired.json_report,
        fastp_paired.html_report
      ]),
      report_title = "PRO-seq Pipeline QC Report"
  }

  output {
    Array[File] fastp_r1_trimmed = fastp_paired.r1_trimmed
    Array[File] fastp_r2_trimmed = fastp_paired.r2_trimmed
    Array[File] fastp_html = fastp_paired.html_report
    Array[File] fastp_json = fastp_paired.json_report
    Array[File?] rrna_unaligned_r1 = deplete_rrna.unaligned_r1
    Array[File?] rrna_unaligned_r2 = deplete_rrna.unaligned_r2
    Array[File] experimental_bam_prededup = align_experimental.sorted_bam
    Array[File] experimental_bam_dedup = dedup_experimental.deduped_bam
    Array[File] experimental_bai_dedup = dedup_experimental.deduped_bai
    Array[File] experimental_dedup_log = dedup_experimental.log
    Array[File] spikein_bam_prededup = keep_spikein_only.filtered_bam
    Array[File] spikein_bam_dedup = dedup_spikein.deduped_bam
    Array[File] spikein_bai_dedup = dedup_spikein.deduped_bai
    Array[File] spikein_dedup_log = dedup_spikein.log
    Array[File] bigwig_forward = bigwig_fwd.coverage_file
    Array[File] bigwig_reverse = bigwig_rev.coverage_file
    File multiqc_report = run_multiqc.html_report
    File multiqc_data = run_multiqc.data_dir
  }
}
