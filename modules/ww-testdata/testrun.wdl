version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/imputation-update/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow testdata_example {
  # Pull down reference genome and index files for chr1
  call ww_testdata.download_ref_data { input:
      chromo = "chr1",
      version = "hg38",
      region = "1-10000000"
  }

  call ww_testdata.download_fastq_data { }

  call ww_testdata.interleave_fastq { input:
    r1_fq = download_fastq_data.r1_fastq,
    r2_fq = download_fastq_data.r2_fastq
  }

  call ww_testdata.download_cram_data { input:
    ref_fasta = download_ref_data.fasta
  }

  call ww_testdata.download_bam_data { }

  call ww_testdata.download_ichor_data { }

  call ww_testdata.download_dbsnp_vcf { input:
    region = "NC_000001.11:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_known_indels_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_gnomad_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_annotsv_vcf { }

  call ww_testdata.generate_pasilla_counts { }

  call ww_testdata.download_test_transcriptome { }

  call ww_testdata.create_clean_amplicon_reference { input:
    input_fasta = download_ref_data.fasta,
    region = "chr1:5000000-5001000",
    output_name = "chr1_test_amplicon",
    replace_n_with = "A"
  }

  call ww_testdata.create_gdc_manifest { }

  call ww_testdata.download_shapemapper_data { }

  call ww_testdata.download_test_cellranger_ref { }

  call ww_testdata.download_diamond_data { }

  call ww_testdata.download_glimpse2_genetic_map { }

  call ww_testdata.download_glimpse2_reference_panel { }

  call ww_testdata.download_glimpse2_test_gl_vcf { }

  call ww_testdata.download_glimpse2_truth_vcf { }

  call validate_outputs { input:
    ref_fasta = download_ref_data.fasta,
    ref_fasta_index = download_ref_data.fasta_index,
    ref_dict = download_ref_data.dict,
    ref_gtf = download_ref_data.gtf,
    ref_bed = download_ref_data.bed,
    r1_fastq = download_fastq_data.r1_fastq,
    r2_fastq = download_fastq_data.r2_fastq,
    inter_fastq = interleave_fastq.inter_fastq,
    cram = download_cram_data.cram,
    crai = download_cram_data.crai,
    bam = download_bam_data.bam,
    bai = download_bam_data.bai,
    ichor_gc_wig = download_ichor_data.wig_gc,
    ichor_map_wig = download_ichor_data.wig_map,
    ichor_centromeres = download_ichor_data.centromeres,
    ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds,
    dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
    dbsnp_vcf_index = download_dbsnp_vcf.dbsnp_vcf_index,
    known_indels_vcf = download_known_indels_vcf.known_indels_vcf,
    known_indels_vcf_index = download_known_indels_vcf.known_indels_vcf_index,
    gnomad_vcf = download_gnomad_vcf.gnomad_vcf,
    gnomad_vcf_index = download_gnomad_vcf.gnomad_vcf_index,
    annotsv_test_vcf = download_annotsv_vcf.test_vcf,
    pasilla_counts = generate_pasilla_counts.individual_count_files,
    pasilla_gene_info = generate_pasilla_counts.gene_info,
    transcriptome_fasta = download_test_transcriptome.transcriptome_fasta,
    clean_amplicon_fasta = create_clean_amplicon_reference.clean_fasta,
    clean_amplicon_fasta_index = create_clean_amplicon_reference.clean_fasta_index,
    clean_amplicon_dict = create_clean_amplicon_reference.clean_dict,
    gdc_manifest = create_gdc_manifest.manifest,
    shapemapper_target_fa = download_shapemapper_data.target_fa,
    shapemapper_modified_r1 = download_shapemapper_data.modified_r1,
    shapemapper_modified_r2 = download_shapemapper_data.modified_r2,
    shapemapper_untreated_r1 = download_shapemapper_data.untreated_r1,
    shapemapper_untreated_r2 = download_shapemapper_data.untreated_r2,
    cellranger_ref_tar = download_test_cellranger_ref.ref_tar,
    diamond_reference = download_diamond_data.reference,
    diamond_query = download_diamond_data.query,
    glimpse2_genetic_map = download_glimpse2_genetic_map.genetic_map,
    glimpse2_reference_vcf = download_glimpse2_reference_panel.reference_vcf,
    glimpse2_reference_vcf_index = download_glimpse2_reference_panel.reference_vcf_index,
    glimpse2_sites_vcf = download_glimpse2_reference_panel.sites_vcf,
    glimpse2_sites_vcf_index = download_glimpse2_reference_panel.sites_vcf_index,
    glimpse2_gl_vcf = download_glimpse2_test_gl_vcf.gl_vcf,
    glimpse2_gl_vcf_index = download_glimpse2_test_gl_vcf.gl_vcf_index,
    glimpse2_truth_vcf = download_glimpse2_truth_vcf.truth_vcf,
    glimpse2_truth_vcf_index = download_glimpse2_truth_vcf.truth_vcf_index
  }

  output {
    # Outputs from the reference data download
    File ref_fasta = download_ref_data.fasta
    File ref_fasta_index = download_ref_data.fasta_index
    File ref_dict = download_ref_data.dict
    File ref_gtf = download_ref_data.gtf
    File ref_bed = download_ref_data.bed
    # Outputs from the fastq, cram, and bam data downloads
    File r1_fastq = download_fastq_data.r1_fastq
    File r2_fastq = download_fastq_data.r2_fastq
    File cram = download_cram_data.cram
    File crai = download_cram_data.crai
    File bam = download_bam_data.bam
    File bai = download_bam_data.bai
    # Outputs from the ichorCNA data download
    File ichor_gc_wig = download_ichor_data.wig_gc
    File ichor_map_wig = download_ichor_data.wig_map
    File ichor_centromeres = download_ichor_data.centromeres
    File ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds
    # Outputs from VCF downloads
    File dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf
    File dbsnp_vcf_index = download_dbsnp_vcf.dbsnp_vcf_index
    File known_indels_vcf = download_known_indels_vcf.known_indels_vcf
    File known_indels_vcf_index = download_known_indels_vcf.known_indels_vcf_index
    File gnomad_vcf = download_gnomad_vcf.gnomad_vcf
    File gnomad_vcf_index = download_gnomad_vcf.gnomad_vcf_index
    File annotsv_test_vcf = download_annotsv_vcf.test_vcf
    # Outputs from Pasilla DESeq2 count generation
    Array[File] pasilla_counts = generate_pasilla_counts.individual_count_files
    Array[String] pasilla_sample_names = generate_pasilla_counts.sample_names
    Array[String] pasilla_sample_conditions = generate_pasilla_counts.sample_conditions
    File pasilla_gene_info = generate_pasilla_counts.gene_info
    # Output from test transcriptome download
    File transcriptome_fasta = download_test_transcriptome.transcriptome_fasta
    # Outputs from clean amplicon reference creation
    File clean_amplicon_fasta = create_clean_amplicon_reference.clean_fasta
    File clean_amplicon_fasta_index = create_clean_amplicon_reference.clean_fasta_index
    File clean_amplicon_dict = create_clean_amplicon_reference.clean_dict
    # Output from GDC manifest creation
    File gdc_manifest = create_gdc_manifest.manifest
    # Outputs from ShapeMapper data download
    File shapemapper_target_fa = download_shapemapper_data.target_fa
    File shapemapper_modified_r1 = download_shapemapper_data.modified_r1
    File shapemapper_modified_r2 = download_shapemapper_data.modified_r2
    File shapemapper_untreated_r1 = download_shapemapper_data.untreated_r1
    File shapemapper_untreated_r2 = download_shapemapper_data.untreated_r2
    # Output from CellRanger reference download
    File cellranger_ref_tar = download_test_cellranger_ref.ref_tar
    # Outputs from DIAMOND data download
    File diamond_reference = download_diamond_data.reference
    File diamond_query = download_diamond_data.query
    # Outputs from GLIMPSE2 test data downloads
    File glimpse2_genetic_map = download_glimpse2_genetic_map.genetic_map
    File glimpse2_reference_vcf = download_glimpse2_reference_panel.reference_vcf
    File glimpse2_reference_vcf_index = download_glimpse2_reference_panel.reference_vcf_index
    File glimpse2_sites_vcf = download_glimpse2_reference_panel.sites_vcf
    File glimpse2_sites_vcf_index = download_glimpse2_reference_panel.sites_vcf_index
    File glimpse2_gl_vcf = download_glimpse2_test_gl_vcf.gl_vcf
    File glimpse2_gl_vcf_index = download_glimpse2_test_gl_vcf.gl_vcf_index
    File glimpse2_truth_vcf = download_glimpse2_truth_vcf.truth_vcf
    File glimpse2_truth_vcf_index = download_glimpse2_truth_vcf.truth_vcf_index
    # Validation report summarizing all outputs
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validates downloaded test data files to ensure they exist and are non-empty"
    outputs: {
        report: "Validation summary reporting file checks and basic statistics"
    }
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA file to validate"
    ref_fasta_index: "Reference FASTA index file to validate"
    ref_dict: "Reference FASTA dictionary file to validate"
    ref_gtf: "GTF annotation file to validate"
    ref_bed: "BED file to validate"
    r1_fastq: "R1 FASTQ file to validate"
    r2_fastq: "R2 FASTQ file to validate"
    inter_fastq: "Interleaved FASTQ to validate"
    cram: "CRAM file to validate"
    crai: "CRAM index file to validate"
    bam: "BAM file to validate"
    bai: "BAM index file to validate"
    ichor_gc_wig: "ichorCNA GC content file to validate"
    ichor_map_wig: "ichorCNA mapping quality file to validate"
    ichor_centromeres: "ichorCNA centromere locations file to validate"
    ichor_panel_of_norm_rds: "ichorCNA panel of normals file to validate"
    dbsnp_vcf: "dbSNP VCF to validate"
    dbsnp_vcf_index: "dbSNP VCF index to validate"
    known_indels_vcf: "Known indels VCF to validate"
    known_indels_vcf_index: "Known indels VCF index to validate"
    gnomad_vcf: "gnomad VCF to validate"
    gnomad_vcf_index: "gnomad VCF index to validate"
    annotsv_test_vcf: "AnnotSV test VCF file to validate"
    pasilla_counts: "Array of individual count files for each sample from Pasilla dataset to validate"
    pasilla_gene_info: "Pasilla gene annotation information to validate"
    transcriptome_fasta: "Test transcriptome FASTA file to validate"
    clean_amplicon_fasta: "Clean amplicon reference FASTA file to validate"
    clean_amplicon_fasta_index: "Clean amplicon reference FASTA index file to validate"
    clean_amplicon_dict: "Clean amplicon reference dictionary file to validate"
    gdc_manifest: "GDC manifest file to validate"
    shapemapper_target_fa: "ShapeMapper target RNA FASTA file to validate"
    shapemapper_modified_r1: "ShapeMapper modified R1 FASTQ file to validate"
    shapemapper_modified_r2: "ShapeMapper modified R2 FASTQ file to validate"
    shapemapper_untreated_r1: "ShapeMapper untreated R1 FASTQ file to validate"
    shapemapper_untreated_r2: "ShapeMapper untreated R2 FASTQ file to validate"
    cellranger_ref_tar: "CellRanger reference tar.gz file to validate"
    diamond_reference: "DIAMOND E. coli reference proteome FASTA file to validate"
    diamond_query: "DIAMOND E. coli query subset FASTA file to validate"
    glimpse2_genetic_map: "GLIMPSE2 genetic map file to validate"
    glimpse2_reference_vcf: "GLIMPSE2 reference panel BCF file to validate"
    glimpse2_reference_vcf_index: "GLIMPSE2 reference panel BCF index to validate"
    glimpse2_sites_vcf: "GLIMPSE2 sites-only VCF file to validate"
    glimpse2_sites_vcf_index: "GLIMPSE2 sites-only VCF index to validate"
    glimpse2_gl_vcf: "GLIMPSE2 genotype likelihoods VCF file to validate"
    glimpse2_gl_vcf_index: "GLIMPSE2 genotype likelihoods VCF index to validate"
    glimpse2_truth_vcf: "GLIMPSE2 truth VCF file to validate"
    glimpse2_truth_vcf_index: "GLIMPSE2 truth VCF index to validate"
    cpu_cores: "Number of CPU cores to use for validation"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_gtf
    File ref_bed
    File r1_fastq
    File r2_fastq
    File inter_fastq
    File cram
    File crai
    File bam
    File bai
    File ichor_gc_wig
    File ichor_map_wig
    File ichor_centromeres
    File ichor_panel_of_norm_rds
    File dbsnp_vcf
    File dbsnp_vcf_index
    File known_indels_vcf
    File known_indels_vcf_index
    File gnomad_vcf
    File gnomad_vcf_index
    File annotsv_test_vcf
    Array[File] pasilla_counts
    File pasilla_gene_info
    File transcriptome_fasta
    File clean_amplicon_fasta
    File clean_amplicon_fasta_index
    File clean_amplicon_dict
    File gdc_manifest
    File shapemapper_target_fa
    File shapemapper_modified_r1
    File shapemapper_modified_r2
    File shapemapper_untreated_r1
    File shapemapper_untreated_r2
    File cellranger_ref_tar
    File diamond_reference
    File diamond_query
    File glimpse2_genetic_map
    File glimpse2_reference_vcf
    File glimpse2_reference_vcf_index
    File glimpse2_sites_vcf
    File glimpse2_sites_vcf_index
    File glimpse2_gl_vcf
    File glimpse2_gl_vcf_index
    File glimpse2_truth_vcf
    File glimpse2_truth_vcf_index
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -euo pipefail

    # Function to validate a file exists and is non-empty
    validate_file() {
      local file_path="$1"
      local file_label="$2"

      if [[ -f "$file_path" && -s "$file_path" ]]; then
        echo "$file_label: $file_path - PASSED" >> validation_report.txt
        return 0
      else
        echo "$file_label: $file_path - MISSING OR EMPTY" >> validation_report.txt
        return 1
      fi
    }

    echo "=== WILDS Test Data Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Validate all files using the function
    validate_file "~{ref_fasta}" "Reference FASTA" || validation_passed=false
    validate_file "~{ref_fasta_index}" "Reference FASTA index" || validation_passed=false
    validate_file "~{ref_dict}" "Reference FASTA dict" || validation_passed=false
    validate_file "~{ref_gtf}" "GTF file" || validation_passed=false
    validate_file "~{ref_bed}" "BED file" || validation_passed=false
    validate_file "~{r1_fastq}" "R1 FASTQ" || validation_passed=false
    validate_file "~{r2_fastq}" "R2 FASTQ" || validation_passed=false
    validate_file "~{inter_fastq}" "Interleaved FASTQ" || validation_passed=false
    validate_file "~{cram}" "CRAM file" || validation_passed=false
    validate_file "~{crai}" "CRAM index" || validation_passed=false
    validate_file "~{bam}" "BAM file" || validation_passed=false
    validate_file "~{bai}" "BAM index" || validation_passed=false
    validate_file "~{ichor_gc_wig}" "ichorCNA GC WIG" || validation_passed=false
    validate_file "~{ichor_map_wig}" "ichorCNA MAP WIG" || validation_passed=false
    validate_file "~{ichor_centromeres}" "ichorCNA centromeres" || validation_passed=false
    validate_file "~{ichor_panel_of_norm_rds}" "ichorCNA panel of normals" || validation_passed=false
    validate_file "~{dbsnp_vcf}" "dbSNP VCF" || validation_passed=false
    validate_file "~{dbsnp_vcf_index}" "dbSNP VCF index" || validation_passed=false
    validate_file "~{known_indels_vcf}" "Known Indels VCF" || validation_passed=false
    validate_file "~{known_indels_vcf_index}" "Known Indels VCF index" || validation_passed=false
    validate_file "~{gnomad_vcf}" "Gnomad VCF" || validation_passed=false
    validate_file "~{gnomad_vcf_index}" "Gnomad VCF index" || validation_passed=false
    validate_file "~{annotsv_test_vcf}" "AnnotSV test VCF" || validation_passed=false
    validate_file "~{cellranger_ref_tar}" "Cellranger test reference"  || validation_passed=false

    # Validate pasilla count files
    for count_file in ~{sep=' ' pasilla_counts}; do
      validate_file "$count_file" "Pasilla count file" || validation_passed=false
    done

    validate_file "~{pasilla_gene_info}" "Pasilla gene info" || validation_passed=false
    validate_file "~{transcriptome_fasta}" "Test transcriptome FASTA" || validation_passed=false
    validate_file "~{gdc_manifest}" "GDC manifest" || validation_passed=false
    validate_file "~{clean_amplicon_fasta}" "Clean amplicon FASTA" || validation_passed=false
    validate_file "~{clean_amplicon_fasta_index}" "Clean amplicon FASTA index" || validation_passed=false
    validate_file "~{clean_amplicon_dict}" "Clean amplicon FASTA dict" || validation_passed=false
    validate_file "~{shapemapper_target_fa}" "ShapeMapper target FASTA" || validation_passed=false
    validate_file "~{shapemapper_modified_r1}" "ShapeMapper modified R1 FASTQ" || validation_passed=false
    validate_file "~{shapemapper_modified_r2}" "ShapeMapper modified R2 FASTQ" || validation_passed=false
    validate_file "~{shapemapper_untreated_r1}" "ShapeMapper untreated R1 FASTQ" || validation_passed=false
    validate_file "~{shapemapper_untreated_r2}" "ShapeMapper untreated R2 FASTQ" || validation_passed=false
    validate_file "~{diamond_reference}" "DIAMOND reference proteome FASTA" || validation_passed=false
    validate_file "~{diamond_query}" "DIAMOND query subset FASTA" || validation_passed=false
    validate_file "~{glimpse2_genetic_map}" "GLIMPSE2 genetic map" || validation_passed=false
    validate_file "~{glimpse2_reference_vcf}" "GLIMPSE2 reference panel BCF" || validation_passed=false
    validate_file "~{glimpse2_reference_vcf_index}" "GLIMPSE2 reference panel BCF index" || validation_passed=false
    validate_file "~{glimpse2_sites_vcf}" "GLIMPSE2 sites-only VCF" || validation_passed=false
    validate_file "~{glimpse2_sites_vcf_index}" "GLIMPSE2 sites-only VCF index" || validation_passed=false
    validate_file "~{glimpse2_gl_vcf}" "GLIMPSE2 genotype likelihoods VCF" || validation_passed=false
    validate_file "~{glimpse2_gl_vcf_index}" "GLIMPSE2 genotype likelihoods VCF index" || validation_passed=false
    validate_file "~{glimpse2_truth_vcf}" "GLIMPSE2 truth VCF" || validation_passed=false
    validate_file "~{glimpse2_truth_vcf_index}" "GLIMPSE2 truth VCF index" || validation_passed=false

    # Additional check: Verify no N bases in clean amplicon
    echo "" >> validation_report.txt
    echo "=== Clean Amplicon Verification ===" >> validation_report.txt
    n_count=$(grep -v "^>" "~{clean_amplicon_fasta}" | grep -o "N" | wc -l | tr -d '[:space:]' || echo "0")
    if [ "$n_count" -eq 0 ]; then
      echo "N base check: PASSED (0 N bases found)" >> validation_report.txt
    else
      echo "N base check: FAILED ($n_count N bases found)" >> validation_report.txt
      validation_passed=false
    fi

    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 45"
    } >> validation_report.txt

    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED - see details above" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
