## WILDS WDL module for downloading FASTQ data from NCBI's Sequence Read Archive (SRA).
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.
## Supports both public SRA data and controlled-access dbGaP data via NGC authentication.

version 1.0

task fastqdump {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for pulling down fastq data from SRA using prefetch and fasterq-dump. Classifies output files by read length to drop index reads and return biological R1/R2 only. Supports controlled-access dbGaP data via optional NGC repository key file."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl"
    outputs: {
        r1_end: "R1 fastq file downloaded for the sample in question",
        r2_end: "R2 fastq file downloaded for the sample in question (empty file for single-end reads)",
        is_paired_end: "boolean indicating whether the sample used paired-end sequencing"
    }
    topic: "any"
    species: "human,eukaryote,prokaryote,virus"
    operation: "data_retrieval"
    input_sample_required: "none"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "r1_end:dna_sequence:fastq,r2_end:dna_sequence:fastq"
    output_reference: "none"
  }

  parameter_meta {
    sra_id: "SRA ID of the sample to be downloaded via prefetch and fasterq-dump"
    ncpu: "number of cpus to use during download"
    max_reads: "Optional maximum number of reads to download (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data. Must be obtained through an approved dbGaP project."
    docker_image: "Docker image to use for this task"
  }

  input {
    String sra_id
    Int ncpu = 8
    Int? max_reads
    File? ngc_file
    String docker_image = "getwilds/sra-tools:3.1.1"
  }

  command <<<
    set -eo pipefail
    # Prefetch the SRA data (handles both public and dbGaP data)
    prefetch "~{sra_id}" \
      ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
    # Pre-create an explicit temp dir for fasterq-dump. Letting it pick
    # its own path has produced "cannot create this temporary directory"
    # failures on NFS-backed execution dirs even when the dir is writable
    # and a shell `mkdir` of the same name succeeds.
    mkdir -p ./fasterq_tmp
    # Check if paired ended
    numLines=$(fastq-dump -X 1 -Z --split-spot "~{sra_id}" \
      ~{if defined(ngc_file) then "--ngc " + ngc_file else ""} | wc -l)
    if [ "$numLines" -gt 4 ]; then
      echo true > paired_file
      fasterq-dump "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --temp ./fasterq_tmp \
        --split-files \
        --include-technical \
        ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
    else
      echo false > paired_file
      fasterq-dump "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --temp ./fasterq_tmp \
        ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
      # Rename the file to match the expected output format
      mv "~{sra_id}.fastq" "~{sra_id}_1.fastq"
      # Create an empty placeholder for R2
      touch "~{sra_id}_2.fastq"
    fi
    # Classify the fasterq-dump output files by read length to figure
    # out which are biological reads vs technical (index) reads. Needed
    # because submissions vary: 10x can produce {R1=28bp, R2=98bp, I1=8bp}
    # in any file-index order, paired bulk + index can produce
    # {76bp, 76bp, 8bp}, and so on. The prior heuristic (always treat _3
    # as R2) silently mis-assigned R1 for 10x submissions whose file order
    # was {R1, R2, I1}, producing an 8bp "R1" that broke Cell Ranger
    # chemistry detection.
    #
    # Rule: any file with median read length <= INDEX_MAX_LEN is treated
    # as an index read and dropped. Of the remaining (biological) files,
    # the one with the longest median is R2; the next is R1. Files of
    # equal median length break ties by original file index, which gives
    # the conventional R1=_1, R2=_2 mapping for plain paired-end data.
    if [ "$numLines" -gt 4 ]; then
      INDEX_MAX_LEN=15
      echo "=== fasterq-dump output classification ==="
      classification_table=""
      for f in "~{sra_id}"_*.fastq; do
        [ -e "$f" ] || continue
        # Median is overkill; first read's length is representative
        # because Illumina runs produce fixed-length reads per file.
        rlen=$(awk 'NR==2 {print length($0); exit}' "$f")
        classification_table="${classification_table}${f}	${rlen}\n"
      done
      printf "file\tread_length\n"
      printf "$classification_table"

      # Split into biological vs index based on length.
      bio_files=()
      bio_lens=()
      for f in "~{sra_id}"_*.fastq; do
        [ -e "$f" ] || continue
        rlen=$(awk 'NR==2 {print length($0); exit}' "$f")
        if [ "$rlen" -le "$INDEX_MAX_LEN" ]; then
          echo "  dropping $f (length $rlen <= $INDEX_MAX_LEN, treated as index read)"
          rm -f "$f"
        else
          bio_files+=("$f")
          bio_lens+=("$rlen")
        fi
      done

      n_bio=${#bio_files[@]}
      if [ "$n_bio" -lt 1 ] || [ "$n_bio" -gt 2 ]; then
        echo "ERROR: expected 1 or 2 biological-read files after dropping indexes, got $n_bio" >&2
        exit 1
      fi

      if [ "$n_bio" -eq 1 ]; then
        # Submission had biological + technical but only one biological
        # read. Treat as single-end: promote the surviving file to R1.
        mv "${bio_files[0]}" "~{sra_id}_1.fastq"
        # Recreate the empty R2 placeholder.
        : > "~{sra_id}_2.fastq"
        # Override the earlier paired_file value since the biological
        # payload is actually single-end.
        echo false > paired_file
        echo "  assigned R1=${bio_files[0]} (single biological read)"
      else
        # Two biological files. Longer-median = R2 (cDNA). Equal lengths
        # tie-break by original index, preserving R1=_1, R2=_2 for plain
        # paired-end submissions.
        if [ "${bio_lens[1]}" -gt "${bio_lens[0]}" ]; then
          r1_src="${bio_files[0]}"
          r2_src="${bio_files[1]}"
        else
          r1_src="${bio_files[1]}"
          r2_src="${bio_files[0]}"
        fi
        # Move via temp names so we don't clobber a source mid-rename.
        mv "$r1_src" "~{sra_id}_r1.tmp.fastq"
        mv "$r2_src" "~{sra_id}_r2.tmp.fastq"
        # Remove any other leftover files (shouldn't be any, but be safe).
        for f in "~{sra_id}"_*.fastq; do
          case "$f" in
            "~{sra_id}_r1.tmp.fastq"|"~{sra_id}_r2.tmp.fastq") ;;
            *) [ -e "$f" ] && rm -f "$f" ;;
          esac
        done
        mv "~{sra_id}_r1.tmp.fastq" "~{sra_id}_1.fastq"
        mv "~{sra_id}_r2.tmp.fastq" "~{sra_id}_2.fastq"
        echo "  assigned R1=$r1_src, R2=$r2_src"
      fi
      echo "=== end classification ==="
    fi
    # Truncate to max_reads if specified (fasterq-dump has no built-in read limit)
    MAX_READS="~{select_first([max_reads, 0])}"
    if [ "$MAX_READS" -gt 0 ]; then
      max_lines=$((MAX_READS * 4))
      for f in "~{sra_id}_1.fastq" "~{sra_id}_2.fastq"; do
        if [ -s "$f" ]; then
          head -n $max_lines "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
        fi
      done
    fi
    # Gzip the output files (fasterq-dump does not support --gzip natively)
    gzip "~{sra_id}_1.fastq" "~{sra_id}_2.fastq"
  >>>

  output {
    File r1_end = "~{sra_id}_1.fastq.gz"
    File r2_end = "~{sra_id}_2.fastq.gz"
    Boolean is_paired_end = read_boolean("paired_file")
  }

  runtime {
    memory: 2 * ncpu + " GB"
    docker: docker_image
    cpu: ncpu
  }
}
