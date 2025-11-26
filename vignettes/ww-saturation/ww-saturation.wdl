version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as bwa_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-saturation/modules/ww-gatk/ww-gatk.wdl" as gatk_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-saturation/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks

struct SaturationSample {
    String name
    File reads
    File? mates
}

workflow saturation_mutagenesis {
    meta {
        author: "Taylor Firman"
        email: "tfirman@fredhutch.org"
        description: "WDL workflow to analyze saturation mutagenesis data using BWA alignment and GATK AnalyzeSaturationMutagenesis"
        url: "https://github.com/getwilds/wilds-wdl-library/vignettes/ww-saturation"
        outputs: {
            variant_counts: "Array of variant count tables for each sample",
            aa_counts: "Array of amino acid count tables for each sample",
            aa_fractions: "Array of amino acid fraction tables for each sample",
            codon_counts: "Array of codon count tables for each sample",
            codon_fractions: "Array of codon fraction tables for each sample",
            cov_length_counts: "Array of coverage length count tables for each sample",
            read_counts: "Array of read count tables for each sample",
            ref_coverage: "Array of reference coverage tables for each sample"
        }
    }

    parameter_meta {
        samples: "List of SaturationSample objects, each containing the sample name, forward reads FASTQ, and optionally reverse reads FASTQs"
        reference_fasta: "Reference genome FASTA file"
        reference_fasta_index: "Index for reference genome FASTA file"
        reference_dict: "Reference genome sequence dictionary"
        orf_range: "Open reading frame range to analyze (e.g., '1-100')"
        cpu_cores: "Number of CPUs to use for BWA alignment and GATK processing"
        memory_gb: "Memory allocation in GB"
    }

    input {
        Array[SaturationSample] samples
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        String orf_range
        Int cpu_cores = 4
        Int memory_gb = 16
    }

    # Build BWA index once for all samples
    call bwa_tasks.bwa_index {
        input:
            reference_fasta = reference_fasta,
            cpu_cores = cpu_cores,
            memory_gb = memory_gb
    }

    # Scatter over the samples
    scatter (sample in samples) {
        # Step 1: Align Nextera sequencing reads to reference using BWA
        call bwa_tasks.bwa_mem {
            input:
                bwa_genome_tar = bwa_index.bwa_index_tar,
                reference_fasta = reference_fasta,
                reads = sample.reads,
                name = sample.name,
                mates = sample.mates,
                paired_end = defined(sample.mates),
                cpu_cores = cpu_cores,
                memory_gb = memory_gb
        }

        # Step 2: Sort BAM by queryname (required for AnalyzeSaturationMutagenesis)
        call samtools_tasks.sort_bam {
            input:
                input_bam = bwa_mem.sorted_bam,
                base_file_name = sample.name,
                cpu_cores = cpu_cores,
                memory_gb = memory_gb
        }

        # Step 3: Analyze Saturation Mutagenesis
        call gatk_tasks.analyze_saturation_mutagenesis {
            input:
                bam = sort_bam.sorted_bam,
                bam_index = sort_bam.sorted_bai,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                orf_range = orf_range,
                base_file_name = sample.name,
                memory_gb = memory_gb,
                cpu_cores = cpu_cores
        }
    }

    output {
        Array[File] variant_counts = analyze_saturation_mutagenesis.variant_counts
        Array[File] aa_counts = analyze_saturation_mutagenesis.aa_counts
        Array[File] aa_fractions = analyze_saturation_mutagenesis.aa_fractions
        Array[File] codon_counts = analyze_saturation_mutagenesis.codon_counts
        Array[File] codon_fractions = analyze_saturation_mutagenesis.codon_fractions
        Array[File] cov_length_counts = analyze_saturation_mutagenesis.cov_length_counts
        Array[File] read_counts = analyze_saturation_mutagenesis.read_counts
        Array[File] ref_coverage = analyze_saturation_mutagenesis.ref_coverage
    }
}
