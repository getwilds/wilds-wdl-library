version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/sra-star-memory/modules/ww-star/ww-star.wdl" as ww_star

struct StarSample {
    String name
    File r1
    File r2
}

workflow star_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_fastq_data { }

  call ww_star.build_index { input:
      reference_fasta = download_ref_data.fasta,
      reference_gtf = download_ref_data.gtf,
      sjdb_overhang = 100,
      genome_sa_index_nbases = 14,
      memory_gb = 8,
      cpu_cores = 2
  }

  # Create test samples array
  Array[StarSample] final_samples = [
    {
      "name": "demo_sample",
      "r1": download_fastq_data.r1_fastq,
      "r2": download_fastq_data.r2_fastq
    }
  ]

  scatter (sample in final_samples) {
    call ww_star.align_two_pass { input:
        star_genome_tar = build_index.star_index_tar,
        r1 = sample.r1,
        r2 = sample.r2,
        name = sample.name,
        sjdb_overhang = 100,
        memory_gb = 8,
        cpu_cores = 2
    }
  }

  call validate_outputs { input:
    bam_files = align_two_pass.bam,
    bai_files = align_two_pass.bai,
    gene_count_files = align_two_pass.gene_counts
  }

  output {
    Array[File] star_bam = align_two_pass.bam
    Array[File] star_bai = align_two_pass.bai
    Array[File] star_gene_counts = align_two_pass.gene_counts
    Array[File] star_log_final = align_two_pass.log_final
    Array[File] star_log_progress = align_two_pass.log_progress
    Array[File] star_log = align_two_pass.log
    Array[File] star_sj = align_two_pass.sj_out
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected STAR output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks and alignment statistics"
    }
  }

  parameter_meta {
    bam_files: "Array of BAM files to validate"
    bai_files: "Array of BAM index files to validate"
    gene_count_files: "Array of gene count files to validate"
  }

  input {
    Array[File] bam_files
    Array[File] bai_files
    Array[File] gene_count_files
  }

  command <<<
    set -eo pipefail
    
    echo "=== STAR Alignment Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt
    
    # Arrays for bash processing
    bam_files=~{sep=" " bam_files}
    bai_files=~{sep=" " bai_files}
    gene_count_files=~{sep=" " gene_count_files}
    
    validation_passed=true
    total_mapped_reads=0
    
    # Check each sample
    for i in "${!bam_files[@]}"; do
      bam_file="${bam_files[$i]}"
      bai_file="${bai_files[$i]}"
      gene_count_file="${gene_count_files[$i]}"
      
      echo "--- Sample: $bam_file ---" >> validation_report.txt
      
      # Check BAM file exists and is not empty
      if [[ -f "$bam_file" && -s "$bam_file" ]]; then
        bam_size=$(stat -c%s "$bam_file")
        echo "BAM file: $bam_file (${bam_size} bytes)" >> validation_report.txt
        
        # Try to get alignment stats from samtools if available
        if command -v samtools &> /dev/null; then
          mapped_reads=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo "N/A")
          total_reads=$(samtools view -c "$bam_file" 2>/dev/null || echo "N/A")
          echo "  Total reads: $total_reads" >> validation_report.txt
          echo "  Mapped reads: $mapped_reads" >> validation_report.txt
          
          if [[ "$mapped_reads" =~ ^[0-9]+$ ]]; then
            total_mapped_reads=$((total_mapped_reads + mapped_reads))
          fi
        fi
      else
        echo "BAM file: $bam_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
      
      # Check BAI file exists
      if [[ -f "$bai_file" ]]; then
        bai_size=$(stat -c%s "$bai_file")
        echo "BAI file: $bai_file (${bai_size} bytes)" >> validation_report.txt
      else
        echo "BAI file: $bai_file - MISSING" >> validation_report.txt
        validation_passed=false
      fi
      
      # Check gene counts file exists and has content
      if [[ -f "$gene_count_file" && -s "$gene_count_file" ]]; then
        gene_count_size=$(stat -c%s "$gene_count_file")
        gene_count_lines=$(wc -l < "$gene_count_file" 2>/dev/null || echo "N/A")
        echo "Gene counts file: $gene_count_file (${gene_count_size} bytes, ${gene_count_lines} lines)" >> validation_report.txt
      else
        echo "Gene counts file: $gene_count_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
      
      echo "" >> validation_report.txt
    done
    
    # Overall summary
    echo "=== Validation Summary ===" >> validation_report.txt
    echo "Total samples processed: ${#bam_files[@]}" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi
    
    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: "2 GB"
    cpu: 1
  }
}
