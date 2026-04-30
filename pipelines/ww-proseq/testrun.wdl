version 1.0

# Local relative imports for the in-PR pipeline + ww-testdata. Flip the testdata
# import back to the GitHub raw URL at merge time (the new helper tasks land in
# the same PR so they don't exist on main yet).
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/pipelines/ww-proseq/ww-proseq.wdl" as proseq_workflow
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-sra/ww-sra.wdl" as ww_sra

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

workflow proseq_example {
  # Pull a small slice of a real qPRO-seq run (Judd et al. 2020, Drosophila S2 cells,
  # GSE149332). Same lab/protocol as the JAJ256/PROseq_alignment.sh source script —
  # 6 bp UMIs on both 5' and 3' ends, fastp + umi_tools dedup workflow.
  call ww_sra.fastqdump as proseq_rep1 { input:
      sra_id = "SRR11607571",
      ncpu = 2,
      max_reads = 50000
  }

  # Experimental genome: a small slice of dm6 (matches the S2-cell biology).
  call ww_testdata.download_ref_data as download_dm6 { input:
      version = "dm6",
      chromo = "chr2L",
      region = "1-1000000",
      output_name = "dm6_chr2L_1M"
  }

  # Spike-in source: a small slice of hg38. The pipeline expects spike-in contigs
  # in the merged reference to share a configurable name prefix.
  call ww_testdata.download_ref_data as download_hg38 { input:
      version = "hg38",
      chromo = "chr1",
      region = "1-1000000",
      output_name = "hg38_chr1_1M"
  }

  # Build the merged experimental + spike-in reference.
  call ww_testdata.merge_fastas_with_prefix as build_merged { input:
      first_fasta = download_dm6.fasta,
      second_fasta = download_hg38.fasta,
      second_prefix = "hg38",
      output_name = "dm6_hg38_merged"
  }

  # rDNA reference for rRNA depletion.
  call ww_testdata.download_rrna_reference { }

  ProseqReferences references = {
      "experimental_fasta": download_dm6.fasta,
      "spikein_merged_fasta": build_merged.merged_fasta,
      "rdna_fasta": download_rrna_reference.fasta,
      "spikein_chrom_prefix": "hg38"
  }

  ProseqSample sample1 = {
      "name": "proseq_demo",
      "r1": proseq_rep1.r1_end,
      "r2": proseq_rep1.r2_end
  }

  call proseq_workflow.proseq { input:
      samples = [sample1],
      references = references,
      umi_loc = "per_read",
      umi_len = 6,
      mapq_threshold = 10,
      align_cpu = 2,
      align_memory_gb = 8
  }

  output {
    Array[File] fastp_r1_trimmed = proseq.fastp_r1_trimmed
    Array[File] fastp_r2_trimmed = proseq.fastp_r2_trimmed
    Array[File] fastp_html = proseq.fastp_html
    Array[File] fastp_json = proseq.fastp_json
    Array[File] experimental_bam_dedup = proseq.experimental_bam_dedup
    Array[File] experimental_bai_dedup = proseq.experimental_bai_dedup
    Array[File] experimental_dedup_log = proseq.experimental_dedup_log
    Array[File] spikein_bam_dedup = proseq.spikein_bam_dedup
    Array[File] spikein_bai_dedup = proseq.spikein_bai_dedup
    Array[File] spikein_dedup_log = proseq.spikein_dedup_log
    Array[File] bigwig_forward = proseq.bigwig_forward
    Array[File] bigwig_reverse = proseq.bigwig_reverse
    File multiqc_report = proseq.multiqc_report
  }
}
