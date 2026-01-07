## Consensus variant calling workflow for DNA sequencing - Refactored with WILDS modules
## This version imports functionality from existing WILDS WDL modules for maximum reusability
##
## Key Features:
## - Modular design using 10 WILDS WDL modules for maintainability and testing
## - Advanced scatter-gather parallelization with interval-based BAM splitting
## - Comprehensive variant calling: germline (HaplotypeCaller), somatic (Mutect2), and population-based (bcftools)
## - Multi-caller consensus variant analysis with custom R processing
## - Complete structural variant analysis pipeline (Manta, Smoove, Delly with AnnotSV annotation)
## - Tumor fraction estimation for cfDNA analysis using ichorCNA
## - GATK best practices for preprocessing (duplicate marking, base quality recalibration)
##
## Input requirements:
## - CRAM files containing paired-end sequencing data
## - Sample information provided as structs in the input JSON
## - Reference genome files (FASTA, index, dictionary)
## - Variant calling resource files (dbSNP, known indels, gnomAD)
## - ichorCNA reference files (GC content, mappability, panel of normals, centromeres)
## - Annovar annotation configuration (protocols and operations)
##
## Major workflow steps:
## 1. CRAM to FASTQ conversion using samtools
## 2. Read alignment using BWA-MEM with proper read groups
## 3. Preprocessing: duplicate marking, base recalibration, and QC metrics
## 4. Interval-based scatter-gather for parallel variant calling
## 5. Three-way variant calling: HaplotypeCaller, Mutect2, and bcftools mpileup
## 6. Variant annotation using Annovar with customizable protocols
## 7. Consensus variant calling combining all three callers
## 8. Structural variant calling and annotation using three complementary tools
## 9. Copy number analysis and tumor fraction estimation using ichorCNA
##
## Output Files:
## - Variant calls from three different callers (VCF format)
## - Annotated variants with functional and clinical annotations
## - Consensus variant calls combining evidence from all callers
## - Structural variant calls with comprehensive genomic annotations
## - Copy number profiles and tumor fraction estimates
## - Quality control metrics and validation reports

version 1.0

# Import WILDS modules
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annotsv/ww-annotsv.wdl" as annotsv_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annovar/ww-annovar.wdl" as annovar_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bcftools/ww-bcftools.wdl" as bcftools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as bwa_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-delly/ww-delly.wdl" as delly_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gatk/ww-gatk.wdl" as gatk_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ichorcna/ww-ichorcna.wdl" as ichorcna_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-manta/ww-manta.wdl" as manta_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-smoove/ww-smoove.wdl" as smoove_tasks

struct SampleDetails {
    String name
    Array[File] cramfiles
}

workflow ww_leukemia {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Consensus variant calling workflow for human panel/PCR-based targeted DNA sequencing with focus on leukemia analysis - Refactored with WILDS modules"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-leukemia/ww-leukemia.wdl"
    outputs: {
        haplotype_vcf: "Array of variant calls from GATK HaplotypeCaller",
        mpileup_vcf: "Array of variant calls from samtools/bcftools mpileup",
        mutect_vcf: "Array of variant calls from GATK Mutect2 (tumor-only mode)",
        mutect_vcf_index: "Array of index files for Mutect2 VCF files",
        mutect_stats: "Array of Mutect2 statistics files summarizing variant calls",
        mutect_annotated_vcf: "Array of Mutect2 VCF files annotated with Annovar",
        mutect_annotated_table: "Array of Mutect2 variants in tabular format with annotations",
        haplotype_annotated_vcf: "Array of GATK HaplotypeCaller VCF files annotated with Annovar",
        haplotype_annotated_table: "Array of GATK HaplotypeCaller variants in tabular format with annotations",
        mpileup_annotated_vcf: "Array of samtools/bcftools VCF files annotated with Annovar",
        mpileup_annotated_table: "Array of samtools/bcftools variants in tabular format with annotations",
        gatk_wgs_metrics: "Array of hybrid selection metrics from GATK CollectWgsMetrics",
        consensus_variants: "Array of consensus variant calls combining results from all three variant callers",
        manta_sv_vcf: "Array of structural variant calls from Manta",
        manta_sv_vcf_index: "Array of index files for the Manta structural variant VCF files",
        manta_sv_annotated_tsv: "Array of Manta SV calls annotated with Annotsv",
        smoove_sv_vcf: "Array of structural variant calls from Smoove",
        smoove_sv_vcf_index: "Array of index files for the Smoove structural variant VCF files",
        smoove_sv_annotated_tsv: "Array of Smoove SV calls annotated with Annotsv",
        delly_sv_bcf: "Array of structural variant calls from Delly",
        delly_sv_bcf_index: "Array of index files for the Delly structural variant BCF files",
        delly_sv_annotated_tsv: "Array of Delly SV calls annotated with Annotsv",
        ichorcna_params: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions",
        ichorcna_seg: "Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male",
        ichorcna_genomewide_pdf: "Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution",
        ichorcna_allgenomewide_pdf: "Combined PDF of all solutions",
        ichorcna_correct_pdf:  "Genome wide correction comparisons",
        ichorcna_rdata: "Saved R image after ichorCNA has finished. Results for all solutions will be included",
        ichorcna_wig: "WIG file created from binned read count data within input BED files"
    }
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file (.fai) for the reference genome FASTA"
    ref_dict: "Sequence dictionary (.dict) for the reference genome"
    dbsnp_vcf: "dbSNP VCF file for base quality score recalibration and variant annotation"
    af_only_gnomad: "gnomAD allele frequency VCF file for Mutect2 germline resource"
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    panel_of_norm_rds: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    ichorcna_chromosomes: "Array of chromosomes for ichorCNA read counting (e.g., ['chr1', 'chr2', ...]). Default: all autosomes + sex chromosomes"
    ichorcna_chrs_string: "R-style string of chromosomes for ichorCNA analysis (e.g., 'c(1:22, \"X\", \"Y\")'). Default: all autosomes + sex chromosomes"
    known_indels_sites_vcfs: "Array of VCF files containing known indel sites for base quality score recalibration"
    samples: "Array of sample information structs containing sample names, molecular IDs, and CRAM files"
    ref_name: "Reference genome build name (e.g., 'hg38', 'hg19') used for Annovar annotation"
    annovar_protocols: "Comma-separated list of Annovar annotation protocols to apply"
    annovar_operation: "Comma-separated list of Annovar operations corresponding to the protocols"
    scatter_count: "Number of intervals to scatter across for parallel processing"
    high_intensity_cpus: "Number of CPU cores for high-intensity tasks (bwa_mem, markdup_recal_metrics, mpileup_call, SV callers). Default: 8 for production, use 2 for testing"
    high_intensity_memory_gb: "Memory allocation in GB for high-intensity tasks. Default: 16 for production, use 4-8 for testing"
    standard_cpus: "Number of CPU cores for standard-intensity tasks (bwa_index, crams_to_fastq, split_intervals, annovar, annotsv, ichorcna). Default: 4 for production, use 2 for testing"
    standard_memory_gb: "Memory allocation in GB for standard-intensity tasks. Default: 8 for production, use 4 for testing"
    skip_annotations: "Skip variant and SV annotation tasks (annovar and annotsv). Set to true for CI/CD testing to reduce disk usage. Default: false (run all annotations)"
  }

  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp_vcf
    File af_only_gnomad
    File wig_gc
    File wig_map
    File panel_of_norm_rds
    File centromeres
    Array[String] ichorcna_chromosomes
    String ichorcna_chrs_string
    Array[File] known_indels_sites_vcfs
    Array[SampleDetails] samples
    String ref_name
    String annovar_protocols
    String annovar_operation
    Int scatter_count = 32
    Int high_intensity_cpus = 8
    Int high_intensity_memory_gb = 16
    Int standard_cpus = 4
    Int standard_memory_gb = 8
    Boolean skip_annotations = false
  }

  # Use GATK module for splitting intervals
  call gatk_tasks.split_intervals { input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = ref_dict,
      scatter_count = scatter_count,
      cpu_cores = standard_cpus,
      memory_gb = standard_memory_gb
  }

  # Use BWA module for indexing
  call bwa_tasks.bwa_index { input:
      reference_fasta = ref_fasta,
      cpu_cores = standard_cpus,
      memory_gb = standard_memory_gb
  }

  scatter (sample in samples){
    String base_file_name = sample.name + "." + ref_name

    # Use Samtools module for CRAM to FASTQ conversion
    call samtools_tasks.crams_to_fastq { input:
        cram_files = sample.cramfiles,
        ref = ref_fasta,
        name = sample.name,
        cpu_cores = standard_cpus,
        memory_gb = standard_memory_gb
    }

    # Use BWA module for alignment
    call bwa_tasks.bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = ref_fasta,
        reads = crams_to_fastq.r1_fastq,
        mates = crams_to_fastq.r2_fastq,
        name = base_file_name,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Use GATK module for mark duplicates and base recalibration (combined)
    call gatk_tasks.markdup_recal_metrics { input:
        bam = bwa_mem.sorted_bam,
        bam_index = bwa_mem.sorted_bai,
        dbsnp_vcf = dbsnp_vcf,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        reference_dict = ref_dict,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        base_file_name = base_file_name,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Run HaplotypeCaller with internal parallelization
    call gatk_tasks.haplotype_caller_parallel { input:
        bam = markdup_recal_metrics.recalibrated_bam,
        bam_index = markdup_recal_metrics.recalibrated_bai,
        intervals = split_intervals.interval_files,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        reference_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        base_file_name = base_file_name,
        memory_gb = 2*scatter_count,
        cpu_cores = scatter_count
    }

    # Run Mutect2 with internal parallelization
    call gatk_tasks.mutect2_parallel { input:
        bam = markdup_recal_metrics.recalibrated_bam,
        bam_index = markdup_recal_metrics.recalibrated_bai,
        intervals = split_intervals.interval_files,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        reference_dict = ref_dict,
        gnomad_vcf = af_only_gnomad,
        base_file_name = base_file_name,
        memory_gb = 2*scatter_count,
        cpu_cores = scatter_count
    }

    # Use bcftools module for variant calling
    call bcftools_tasks.mpileup_call { input:
        bam_file = markdup_recal_metrics.recalibrated_bam,
        bam_index = markdup_recal_metrics.recalibrated_bai,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Use Annovar module for variant annotation (skippable for CI/CD testing)
    if (!skip_annotations) {
      call annovar_tasks.annovar_annotate as annotateSAM { input:
          vcf_to_annotate = mpileup_call.mpileup_vcf,
          ref_name = ref_name,
          annovar_operation = annovar_operation,
          annovar_protocols = annovar_protocols,
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }

      call annovar_tasks.annovar_annotate as annotateMutect { input:
          vcf_to_annotate = mutect2_parallel.vcf,
          ref_name = ref_name,
          annovar_operation = annovar_operation,
          annovar_protocols = annovar_protocols,
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }

      call annovar_tasks.annovar_annotate as annotateHaplotype { input:
          vcf_to_annotate = haplotype_caller_parallel.vcf,
          ref_name = ref_name,
          annovar_operation = annovar_operation,
          annovar_protocols = annovar_protocols,
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }

      # Keep custom consensus processing task
      call consensus_processing { input:
          gatk_vars = annotateHaplotype.annotated_table,
          sam_vars = annotateSAM.annotated_table,
          mutect_vars = annotateMutect.annotated_table,
          base_file_name = base_file_name
      }
    }

    # Use Manta module for structural variants
    call manta_tasks.manta_call { input:
        aligned_bam = markdup_recal_metrics.recalibrated_bam,
        aligned_bam_index = markdup_recal_metrics.recalibrated_bai,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        sample_name = base_file_name,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Use Smoove module for structural variants
    call smoove_tasks.smoove_call { input:
        aligned_bam = markdup_recal_metrics.recalibrated_bam,
        aligned_bam_index = markdup_recal_metrics.recalibrated_bai,
        sample_name = base_file_name,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Use Delly module for structural variants
    call delly_tasks.delly_call { input:
        aligned_bam = markdup_recal_metrics.recalibrated_bam,
        aligned_bam_index = markdup_recal_metrics.recalibrated_bai,
        reference_fasta = ref_fasta,
        reference_fasta_index = ref_fasta_index,
        cpu_cores = high_intensity_cpus,
        memory_gb = high_intensity_memory_gb
    }

    # Use AnnotSV module for SV annotation (skippable for CI/CD testing)
    if (!skip_annotations) {
      call annotsv_tasks.annotsv_annotate as annotateManta { input:
          raw_vcf = manta_call.vcf,
          genome_build = if ref_name == "hg38" then "GRCh38" else "GRCh37",
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }

      call annotsv_tasks.annotsv_annotate as annotateSmoove { input:
          raw_vcf = smoove_call.vcf,
          genome_build = if ref_name == "hg38" then "GRCh38" else "GRCh37",
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }

      call annotsv_tasks.annotsv_annotate as annotateDelly { input:
          raw_vcf = delly_call.vcf,
          genome_build = if ref_name == "hg38" then "GRCh38" else "GRCh37",
          cpu_cores = standard_cpus,
          memory_gb = standard_memory_gb
      }
    }

    # Use ichorCNA module for tumor fraction analysis
    call ichorcna_tasks.readcounter_wig { input:
        bam_file = markdup_recal_metrics.recalibrated_bam,
        bam_index = markdup_recal_metrics.recalibrated_bai,
        sample_name = base_file_name,
        chromosomes = ichorcna_chromosomes,
        window_size = 500000,
        cpus = standard_cpus,
        memory_gb = standard_memory_gb
    }

    call ichorcna_tasks.ichorcna_call { input:
      wig_tumor = readcounter_wig.wig_file,
      wig_gc = wig_gc,
      wig_map = wig_map,
      panel_of_norm_rds = panel_of_norm_rds,
      centromeres = centromeres,
      name = base_file_name,
      sex = "male", # Defaulting to male for now, can be parameterized later
      genome = ref_name,
      genome_style = "UCSC", # Defaulting to UCSC style for now
      chrs = ichorcna_chrs_string,
      cpus = standard_cpus,
      memory_gb = standard_memory_gb
    }
  } # End scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] haplotype_vcf = haplotype_caller_parallel.vcf
    Array[File] mpileup_vcf = mpileup_call.mpileup_vcf
    Array[File] mutect_vcf = mutect2_parallel.vcf
    Array[File] mutect_vcf_index = mutect2_parallel.vcf_index
    Array[File] mutect_stats = mutect2_parallel.stats_file
    Array[File] mutect_annotated_vcf = select_all(annotateMutect.annotated_vcf)
    Array[File] mutect_annotated_table = select_all(annotateMutect.annotated_table)
    Array[File] haplotype_annotated_vcf = select_all(annotateHaplotype.annotated_vcf)
    Array[File] haplotype_annotated_table = select_all(annotateHaplotype.annotated_table)
    Array[File] mpileup_annotated_vcf = select_all(annotateSAM.annotated_vcf)
    Array[File] mpileup_annotated_table = select_all(annotateSAM.annotated_table)
    Array[File] gatk_wgs_metrics = markdup_recal_metrics.wgs_metrics
    Array[File] consensus_variants = select_all(consensus_processing.consensus_tsv)
    Array[File] manta_sv_vcf = manta_call.vcf
    Array[File] manta_sv_vcf_index = manta_call.vcf_index
    Array[File] manta_sv_annotated_tsv = select_all(annotateManta.annotated_tsv)
    Array[File] smoove_sv_vcf = smoove_call.vcf
    Array[File] smoove_sv_vcf_index = smoove_call.vcf_index
    Array[File] smoove_sv_annotated_tsv = select_all(annotateSmoove.annotated_tsv)
    Array[File] delly_sv_bcf = delly_call.vcf
    Array[File] delly_sv_bcf_index = delly_call.vcf_index
    Array[File] delly_sv_annotated_tsv = select_all(annotateDelly.annotated_tsv)
    Array[File] ichorcna_params = ichorcna_call.params
    Array[File] ichorcna_seg = ichorcna_call.seg
    Array[File] ichorcna_genomewide_pdf = ichorcna_call.genomewide_pdf
    Array[File] ichorcna_allgenomewide_pdf = ichorcna_call.allgenomewide_pdf
    Array[File] ichorcna_correct_pdf = ichorcna_call.correct_pdf
    Array[File] ichorcna_rdata = ichorcna_call.rdata
    Array[File] ichorcna_wig = readcounter_wig.wig_file
  }
}

#### CUSTOM TASK DEFINITIONS (Tasks not available in WILDS modules)

# Keep custom consensus processing task (specialized R script)
task consensus_processing {
  meta {
    description: "Generate consensus variant calls by combining results from multiple variant callers"
    outputs: {
        consensus_tsv: "Tab-separated file containing consensus variant calls with evidence from all callers"
    }
  }

  parameter_meta {
    gatk_vars: "Annotated variant table from GATK HaplotypeCaller"
    sam_vars: "Annotated variant table from samtools/bcftools"
    mutect_vars: "Annotated variant table from GATK Mutect2"
    base_file_name: "Base name for output files"
  }

  input {
    File gatk_vars
    File sam_vars
    File mutect_vars
    String base_file_name
  }

  command <<<
    set -eo pipefail
    Rscript /consensus-trio-unpaired.R \
      "~{gatk_vars}" "~{sam_vars}" "~{mutect_vars}" "~{base_file_name}"
  >>>

  output {
    File consensus_tsv = "~{base_file_name}.consensus.tsv"
  }

  runtime {
    cpu: 1
    memory: "8 GB"
    docker: "getwilds/consensus:0.1.1"
  }
}
