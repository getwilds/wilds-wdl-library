## WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File counts_bed
}

workflow ichorcna_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for cfDNA tumor fraction estimation via ichorCNA"
    url: "https://github.com/getwilds/ww-ichorcna"
    outputs: {
        params: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions",
        seg: "Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male",
        genomewide_pdf: "Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution",
        allgenomewide_pdf: "Combined PDF of all solutions",
        correct_pdf:  "Genome wide correction comparisons",
        rdata: "Saved R image after ichorCNA has finished. Results for all solutions will be included",
        wig: "WIG file created from binned read count data within input BED files",
    }
  }

  parameter_meta {
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    norm_depth: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    samples: "Array of sample information containing name and tarball of per-chromosome BED files of read counts"
    sex: "User-specified: male or female"
    genome: "Genome build (e.g. hg19)"
    genome_style: "Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC"
    memory_gb: "Memory allocated for each task in the workflow in GB"
    cpus: "Number of CPU cores allocated for each task in the workflow"
  }

  input {
    File wig_gc
    File wig_map
    File norm_depth
    File centromeres
    Array[SampleInfo] samples
    String sex
    String genome
    String genome_style
    Int memory_gb
    Int cpus
  }

  scatter (sample in samples) {
    call ichorcna_call { input:
          wig_gc = wig_gc,
      wig_map = wig_map,
      norm_depth = norm_depth,
      centromeres = centromeres,
      counts_bed = sample.counts_bed,
      name = sample.name,
      sex = sex,
      genome = genome,
      genome_style = genome_style,
      cpus = cpus,
      memory_gb = memory_gb,
    }
  }

    # call validate_outputs { input:
    #     vcf_files = manta_call.vcf,
    #     vcf_index_files = manta_call.vcf_index
    # }

  output {
    Array[File] params = ichorcna_call.params
    Array[File] seg = ichorcna_call.seg
    Array[File] genomewide_pdf = ichorcna_call.genomewide_pdf
    Array[File] allgenomewide_pdf = ichorcna_call.allgenomewide_pdf
    Array[File] correct_pdf = ichorcna_call.correct_pdf
    Array[File] rdata = ichorcna_call.rdata
    Array[File] wig = ichorcna_call.wig
  }
}

task ichorcna_call {
  meta {
    description: "Estimate cfDNA tumor fraction using ichorCNA"
    outputs: {
        params: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions",
        seg: "Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male",
        genomewide_pdf: "Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution",
        allgenomewide_pdf: "Combined PDF of all solutions",
        correct_pdf:  "Genome wide correction comparisons",
        rdata: "Saved R image after ichorCNA has finished. Results for all solutions will be included",
        wig: "WIG file created from binned read count data within input BED files",
    }
  }

  parameter_meta {
    counts_bed: "Tarball of per-chromosome BED files of read counts. Used to make tumor WIG fileq"
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    norm_depth: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    name: "Sample ID"
    sex: "User-specified: male or female"
    genome: "Genome build (e.g. hg19)"
    genome_style: "Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC"
    memory_gb: "Memory allocated for each task in the workflow in GB"
    cpus: "Number of CPU cores allocated for each task in the workflow"
  }

  input {
    File counts_bed
    File wig_gc
    File wig_map
    File norm_depth
    File centromeres
    String name
    String sex
    String genome
    String genome_style
    Int memory_gb
    Int cpus
  }

  command <<<
    set -eo pipefail
    # Extract tarball of BED files and concatenate
    tar -xzf "~{counts_bed}" -O | sort -k 1,1V -k 2,2n | \
    # Create fixed-step WIG format using 500kb windows
    awk -v window=500000 \
    'BEGIN { chr=""; } { if ($1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\n",$1,window,window); chr=$1; } print $4; }' \
    > "~{name}.ichor.tumor.wig" && \
    /usr/local/bin/Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R \
    --id "~{name}" \
    --WIG "~{name}.ichor.tumor.wig" \
    --ploidy "c(2)" \
    --normal "c(0.1,0.5,.85)" \
    --maxCN 3 \
    --gcWig "~{wig_gc}" \
    --mapWig "~{wig_map}" \
    --centromere "~{centromeres}" \
    --normalPanel "~{norm_depth}" \
    --genomeBuild "~{genome}" \
    --sex "~{sex}" \
    --chrs "c(1:22, \"X\", \"Y\")" \
    --fracReadsInChrYForMale 0.0005 \
    --txnE 0.999999 \
    --txnStrength 1000000 \
    --genomeStyle "~{genome_style}" \
    --libdir /usr/local/bin/ichorCNA/ && \
    # Keep only biologically relevant segments by keeping all Y segments if male
    # and only non-Y segments if female
    awk -v G="~{sex}" '$2!~/Y/ || G=="male"' "~{name}.seg.txt" \
    > "~{name}.ichor.segs.txt" && \
    mv "~{name}"/*.pdf .
  >>>

  output {
    File params = "${name}.params.txt"
    File seg = "${name}.ichor.segs.txt"
    File genomewide_pdf = "${name}_genomeWide.pdf"
    File allgenomewide_pdf = "${name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${name}_genomeWideCorrection.pdf"
    File rdata = "${name}.RData"
    File wig = "${name}.ichor.tumor.wig"
  }

  runtime {
    docker: "fredhutch/ichorcna:0.5.0"
    cpu: cpus
    memory: "~{memory_gb} GB"
  }
}

# task validate_outputs {
#   meta {
#     description: "Validate Manta outputs and generate a comprehensive report"
#     outputs: {
#         report: "Validation summary with structural variant calling statistics"
#     }
#   }

#   parameter_meta {
#     vcf_files: "Array of VCF files to validate"
#     vcf_index_files: "Array of VCF index files to validate"
#   }

#   input {
#     Array[File] vcf_files
#     Array[File] vcf_index_files
#   }

#   command <<<
#     set -eo pipefail

#     echo "Manta Structural Variant Calling Validation Report" > validation_report.txt
#     echo "=================================================" >> validation_report.txt
#     echo "Generated on: $(date)" >> validation_report.txt
#     echo "" >> validation_report.txt

#     echo "Sample Summary:" >> validation_report.txt
#     echo "Total samples processed: ~{length(vcf_files)}" >> validation_report.txt
#     echo "" >> validation_report.txt

#     # Validate each sample's outputs
#     vcf_files=(~{sep=" " vcf_files})
#     vcf_index_files=(~{sep=" " vcf_index_files})

#     for i in "${!vcf_files[@]}"; do
#         vcf="${vcf_files[$i]}"
#         vcf_index="${vcf_index_files[$i]}"

#         echo "Sample: $vcf" >> validation_report.txt

#         # Check if VCF file exists and is not empty
#         if [[ -f "$vcf" && -s "$vcf" ]]; then
#             echo "VCF file present and non-empty" >> validation_report.txt
#             variant_count=$(zcat "$vcf" | grep -v '^#' | wc -l)
#             echo "Variants called: $variant_count" >> validation_report.txt
#         else
#             echo "VCF file missing or empty" >> validation_report.txt
#         fi

#         # Check if VCF index exists
#         if [[ -f "$vcf_index" ]]; then
#             echo "VCF index file present" >> validation_report.txt
#         else
#             echo "VCF index file missing" >> validation_report.txt
#         fi

#         echo "" >> validation_report.txt
#     done

#     echo "Validation completed successfully" >> validation_report.txt
#   >>>

#   output {
#     File report = "validation_report.txt"
#   }

#   runtime {
#     docker: "getwilds/manta:1.6.0"
#     memory: "4GB"
#     cpu: 1
#   }
# }
