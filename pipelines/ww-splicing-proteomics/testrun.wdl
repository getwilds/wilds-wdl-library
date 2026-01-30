version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/pipelines/ww-splicing-proteomics/ww-splicing-proteomics.wdl" as splicing_proteomics_workflow

struct SampleInfo {
    String name
    File r1
    File r2
    String group
}

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow splicing_proteomics_example {
  meta {
    description: "Test workflow for the splicing proteomics pipeline using airway RNA-seq data"
  }

  # Download Ensembl chr15 reference data from JCAST test data
  # This provides both the Ensembl GTF (required by JCAST) and genome FASTA
  call ww_testdata.download_jcast_test_data { }

  # Download RNA-seq data from SRA (airway study - dexamethasone treatment)
  # Using 2 treated + 2 untreated samples for differential splicing analysis
  # Using 500K reads per sample to have sufficient coverage for splicing detection
  call ww_sra.fastqdump as untreated1 { input: sra_id = "SRR1039509", ncpu = 2, max_reads = 500000 }
  call ww_sra.fastqdump as untreated2 { input: sra_id = "SRR1039513", ncpu = 2, max_reads = 500000 }
  call ww_sra.fastqdump as treated1 { input: sra_id = "SRR1039508", ncpu = 2, max_reads = 500000 }
  call ww_sra.fastqdump as treated2 { input: sra_id = "SRR1039512", ncpu = 2, max_reads = 500000 }

  # Construct SampleInfo structs from SRA downloads
  # group1 = untreated (control), group2 = treated (dexamethasone)
  SampleInfo sample1 = {
    "name": "untreated_1",
    "r1": untreated1.r1_end,
    "r2": untreated1.r2_end,
    "group": "group1"
  }

  SampleInfo sample2 = {
    "name": "untreated_2",
    "r1": untreated2.r1_end,
    "r2": untreated2.r2_end,
    "group": "group1"
  }

  SampleInfo sample3 = {
    "name": "treated_1",
    "r1": treated1.r1_end,
    "r2": treated1.r2_end,
    "group": "group2"
  }

  SampleInfo sample4 = {
    "name": "treated_2",
    "r1": treated2.r1_end,
    "r2": treated2.r2_end,
    "group": "group2"
  }

  # Construct RefGenome struct using Ensembl chr15 data
  RefGenome reference = {
    "name": "GRCh38_chr15",
    "fasta": download_jcast_test_data.genome_fasta,
    "gtf": download_jcast_test_data.gtf_file
  }

  # Run splicing proteomics pipeline
  call splicing_proteomics_workflow.splicing_proteomics {
    input:
      samples = [sample1, sample2, sample3, sample4],
      reference_genome = reference,
      ensembl_gtf = download_jcast_test_data.gtf_file,
      read_length = 63,  # Airway dataset read length
      read_type = "paired",
      library_type = "fr-unstranded",
      output_prefix = "airway_splicing",
      sjdb_overhang = 62,  # read_length - 1
      genome_sa_index_nbases = 11,  # Reduced for chr15 subset
      star_cpu = 4,
      star_memory_gb = 16,
      rmats_cpu = 4,
      rmats_memory_gb = 8,
      jcast_cpu = 2,
      jcast_memory_gb = 4
  }

  output {
    Array[File] star_bam = splicing_proteomics.star_bam
    Array[File] star_bai = splicing_proteomics.star_bai
    Array[File] star_gene_counts = splicing_proteomics.star_gene_counts
    Array[File] star_log_final = splicing_proteomics.star_log_final
    Array[File] star_sj = splicing_proteomics.star_sj
    File rmats_output = splicing_proteomics.rmats_output
    File jcast_protein_fasta = splicing_proteomics.jcast_protein_fasta
    File jcast_output = splicing_proteomics.jcast_output
  }
}
