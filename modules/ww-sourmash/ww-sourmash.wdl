## WILDS WDL Module: ww-sourmash
## Purpose: Generate MinHash signatures, perform similarity searches, and metagenome analysis
## Documentation: https://sourmash.readthedocs.io/

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct SourmashSample {
    String name
    File fastq_file
}

workflow sourmash_example {
    meta {
        description: "Demonstration workflow for sourmash module functionality including sketch, search, and gather operations"
        outputs: {
            query_signatures: "MinHash signatures generated from query sequences",
            database_signatures: "MinHash signatures generated from database sequences",
            search_results: "Similarity search results between query and database signatures",
            gather_results: "Metagenome gather analysis results identifying constituent genomes",
            validation_report: "Comprehensive validation report for all outputs"
        }
    }

    parameter_meta {
        samples: "Array of sample information with names and FASTQ files"
        database_fastas: "Array of database FASTA files for signature generation"
        ksize: "k-mer size for signature generation (recommended: 21, 31, 51)"
        threshold: "Minimum threshold for reporting search matches (0.0-1.0)"
        molecule_type: "Molecule type: dna, protein, dayhoff, hp"
        scaled: "Use scaled MinHash (recommended: true for most use cases)"
        scale_factor: "Scale factor for scaled MinHash (lower = more precise, higher = more memory efficient)"
        gather_threshold: "Minimum threshold for gather analysis (0.0-1.0)"
    }

    input {
        Array[SourmashSample]? samples
        Array[File]? database_fastas
        String ksize = "31"
        Float threshold = 0.08
        String molecule_type = "dna"
        Boolean scaled = true
        Int scale_factor = 1000
        Float gather_threshold = 0.05
    }

    # Download test data if no samples provided
    if (!defined(samples)) {
        call ww_testdata.download_fastq_data { }
    }

    # Create samples array - either from input or from ww-testdata
    Array[SourmashSample] final_samples = if defined(samples) then select_first([samples]) else [
      {
        "name": "demo_sample1",
        "fastq_file": select_first([download_fastq_data.r1_fastq])
      },
      {
        "name": "demo_sample2",
        "fastq_file": select_first([download_fastq_data.r2_fastq])
      }
    ]

    # Download test database if no database files provided
    if (!defined(database_fastas)) {
        call ww_testdata.download_ref_data { }
    }

    # Use provided fastas or test reference data
    Array[File] final_fastas = select_first([database_fastas, [select_first([download_ref_data.fasta])]])

    # Generate signatures for all sample files
    scatter (sample_to_sketch in final_samples) {
        call sourmash_sketch as sketch_samples { input:
            input_file = sample_to_sketch.fastq_file,
            output_name = sample_to_sketch.name + ".sig",
            ksize = ksize,
            molecule_type = molecule_type,
            scaled = scaled,
            scale_factor = scale_factor
        }
    }

    # Generate signatures for all database files  
    scatter (fasta in final_fastas) {
        call sourmash_sketch as sketch_database { input:
            input_file = fasta,
            output_name = "database_" + basename(basename(fasta, ".fa"),".fasta") + ".sig",
            ksize = ksize,
            molecule_type = molecule_type,
            scaled = scaled,
            scale_factor = scale_factor
        }
    }

    scatter (sample_to_search in sketch_samples.signature_file) {
        scatter (database in sketch_database.signature_file) {
            # Perform search for each sample against each database signature
            call sourmash_search as search_samples { input:
                query_sig = sample_to_search,
                database_sig = database,
                ksize = ksize,
                threshold = threshold
            }
        }

        # Perform gather analysis for each sample against the entire database
        call sourmash_gather as gather_samples { input:
            query_sig = sample_to_search,
            database_sigs = sketch_database.signature_file,
            ksize = ksize,
            threshold = gather_threshold
        }
    }

    # Validate all outputs
    call validate_outputs { input:
        query_signatures = sketch_samples.signature_file,
        database_signatures = sketch_database.signature_file,
        search_results = flatten(search_samples.results_file),
        gather_results = gather_samples.results_file
    }

    output {
        Array[File] query_signatures = sketch_samples.signature_file
        Array[File] database_signatures = sketch_database.signature_file
        Array[Array[File]] search_results = search_samples.results_file
        Array[File] gather_results = gather_samples.results_file
        File validation_report = validate_outputs.report
    }
}

task sourmash_sketch {
    meta {
        description: "Generate MinHash signatures from FASTA/FASTQ files using sourmash"
        outputs: {
            signature_file: "MinHash signature file in JSON format"
        }
    }

    parameter_meta {
        input_file: "Input sequence file (FASTA or FASTQ format, optionally compressed)"
        output_name: "Name for the output signature file"
        ksize: "k-mer size for signature generation"
        molecule_type: "Molecule type: dna, protein, dayhoff, hp, or nucleotide"
        scaled: "Use scaled MinHash (recommended for most applications)"
        scale_factor: "Scale factor for scaled MinHash (lower = more precise, higher = faster)"
        memory_gb: "Memory allocation in GB"
        cpu_cores: "Number of CPU cores to use"
    }

    input {
        File input_file
        String output_name
        String ksize
        String molecule_type
        Boolean scaled
        Int scale_factor
        Int memory_gb = 4
        Int cpu_cores = 2
    }

    command <<<
        set -eo pipefail
        
        # Generate MinHash signature
        sourmash sketch \
            ~{molecule_type} \
            ~{input_file} \
            -p "k=~{ksize},~{if scaled then "scaled=" else "num="}~{scale_factor}" \
            -o ~{output_name}
    >>>

    output {
        File signature_file = output_name
    }

    runtime {
        docker: "getwilds/sourmash:4.8.2"
        memory: "~{memory_gb}GB"
        cpu: cpu_cores
    }
}

task sourmash_search {
    meta {
        description: "Compare MinHash signatures to find sequence similarities using containment analysis"
        outputs: {
            results_file: "CSV file containing similarity search results"
        }
    }

    parameter_meta {
        query_sig: "Query MinHash signature file"
        database_sig: "Database MinHash signature file to search against"
        ksize: "k-mer size to use for search (must match signature k-mer size)"
        threshold: "Minimum similarity threshold for reporting matches (0.0-1.0)"
        memory_gb: "Memory allocation in GB"
        cpu_cores: "Number of CPU cores to use"
    }

    input {
        File query_sig
        File database_sig
        String ksize
        Float threshold
        Int memory_gb = 4
        Int cpu_cores = 1
    }

    command <<<
        set -eo pipefail

        # Perform containment search between signatures
        sourmash search \
            -k ~{ksize} \
            --threshold ~{threshold} \
            --containment \
            ~{query_sig} \
            ~{database_sig} \
            > search_results.csv
    >>>

    output {
        File results_file = "search_results.csv"
    }

    runtime {
        docker: "getwilds/sourmash:4.8.2"
        memory: "~{memory_gb}GB"
        cpu: cpu_cores
    }
}

task sourmash_gather {
    meta {
        description: "Analyze metagenome samples by identifying constituent genomes using sourmash gather"
        outputs: {
            results_file: "CSV file containing gather analysis results with genome matches and abundances"
        }
    }

    parameter_meta {
        query_sig: "Query MinHash signature file (typically from a metagenome sample)"
        database_sigs: "Array of database MinHash signature files representing reference genomes"
        ksize: "k-mer size to use for gather analysis (must match signature k-mer size)"
        threshold: "Minimum abundance threshold for reporting genome matches (0.0-1.0)"
        memory_gb: "Memory allocation in GB"
        cpu_cores: "Number of CPU cores to use"
    }

    input {
        File query_sig
        Array[File] database_sigs
        String ksize
        Float threshold
        Int memory_gb = 8
        Int cpu_cores = 2
    }

    command <<<
        set -eo pipefail

        # Perform gather analysis to identify constituent genomes
        sourmash gather \
            -k ~{ksize} \
            --threshold ~{threshold} \
            ~{query_sig} \
            ~{sep=' ' database_sigs} \
            -o gather_results.csv
    >>>

    output {
        File results_file = "gather_results.csv"
    }

    runtime {
        docker: "getwilds/sourmash:4.8.2"
        memory: "~{memory_gb}GB"
        cpu: cpu_cores
    }
}

task validate_outputs {
    meta {
        description: "Validate sourmash outputs including signatures and analysis results"
        outputs: {
            report: "Validation report confirming output integrity and providing basic statistics"
        }
    }

    parameter_meta {
        query_signatures: "Query MinHash signature files to validate"
        database_signatures: "Database MinHash signature files to validate"
        search_results: "Search results CSV file to validate"
        gather_results: "Gather results CSV file to validate"
        cpu_cores: "Number of CPU cores to use for validation"
        memory_gb: "Memory allocation in GB for the task"
    }

    input {
        Array[File] query_signatures
        Array[File] database_signatures
        Array[File] search_results
        Array[File] gather_results
        Int cpu_cores = 1
        Int memory_gb = 2
    }

    command <<<
        set -euo pipefail

        echo "=== WILDS Sourmash Module Validation Report ===" > validation_report.txt
        echo "Generated on: $(date)" >> validation_report.txt
        echo "" >> validation_report.txt

        validation_passed=true

        # Validate query signatures
        echo "Query Signatures:" >> validation_report.txt
        ~{sep=' ' query_signatures} | tr ' ' '\n' | while read sig_file; do
            if [[ -f "$sig_file" && -s "$sig_file" ]]; then
                # Count signatures in file
                sig_count=$(sourmash sig info "$sig_file" 2>/dev/null | grep -c "signature" || echo "0")
                echo "  $sig_file - PASSED ($sig_count signatures)" >> validation_report.txt
            else
                echo "  $sig_file - MISSING OR EMPTY" >> validation_report.txt
                validation_passed=false
            fi
        done

        # Validate database signatures
        echo "" >> validation_report.txt
        echo "Database Signatures:" >> validation_report.txt
        ~{sep=' ' database_signatures} | tr ' ' '\n' | while read sig_file; do
            if [[ -f "$sig_file" && -s "$sig_file" ]]; then
                sig_count=$(sourmash sig info "$sig_file" 2>/dev/null | grep -c "signature" || echo "0")
                echo "  $sig_file - PASSED ($sig_count signatures)" >> validation_report.txt
            else
                echo "  $sig_file - MISSING OR EMPTY" >> validation_report.txt
                validation_passed=false
            fi
        done

        # Validate search results
        echo "" >> validation_report.txt
        echo "Search Results:" >> validation_report.txt
        ~{sep=' ' search_results} | tr ' ' '\n' | while read search_file; do
            if [[ -f "$search_file" && -s "$search_file" ]]; then
                result_count=$(tail -n +2 "$search_file" | wc -l)
                echo "  $search_file - PASSED ($result_count matches found)" >> validation_report.txt
            else
                echo "  $search_file - MISSING OR EMPTY" >> validation_report.txt
                validation_passed=false
            fi
        done

        # Validate gather results
        echo "" >> validation_report.txt
        echo "Gather Results:" >> validation_report.txt
        ~{sep=' ' gather_results} | tr ' ' '\n' | while read gather_file; do
            if [[ -f "$gather_file" && -s "$gather_file" ]]; then
                result_count=$(tail -n +2 "$gather_file" | wc -l)
                echo "  $gather_file - PASSED ($result_count genome matches found)" >> validation_report.txt
            else
                echo "  $gather_file - MISSING OR EMPTY" >> validation_report.txt
                validation_passed=false
            fi
        done

        # Final validation status
        echo "" >> validation_report.txt
        if [[ "$validation_passed" == "true" ]]; then
            echo "Overall Status: VALIDATION PASSED - All files present and valid" >> validation_report.txt
        else
            echo "Overall Status: VALIDATION FAILED - Some files missing or invalid" >> validation_report.txt
            exit 1
        fi
    >>>

    output {
        File report = "validation_report.txt"
    }

    runtime {
        docker: "getwilds/sourmash:4.8.2"
        memory: "~{memory_gb}GB"
        cpu: cpu_cores
    }
}
