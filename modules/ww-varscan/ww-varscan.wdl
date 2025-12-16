
## WILDS WDL module for VarScan variant calling.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task somatic {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run VarScan somatic variant calling on tumor-normal pair"
    outputs: {
        somatic_snvs_vcf: "VCF file containing somatic SNV calls",
        somatic_indels_vcf: "VCF file containing somatic indel calls"
    }
  }

  parameter_meta {
    sample_name: "Name of the sample (used in output file names)"
    normal_pileup: "Samtools mpileup file for the normal sample"
    tumor_pileup: "Samtools mpileup file for the tumor sample"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    String sample_name
    File normal_pileup
    File tumor_pileup
    Int memory_gb = 16
    Int cpu_cores = 4
  }

  command <<<
    set -euo pipefail

    java -Xmx~{memory_gb - 2}g -jar /usr/local/bin/VarScan.jar somatic \
      "~{normal_pileup}" \
      "~{tumor_pileup}" \
      "~{sample_name}" \
      --output-vcf 1
  >>>

  output {
    File somatic_snvs_vcf = "~{sample_name}.snp.vcf"
    File somatic_indels_vcf = "~{sample_name}.indel.vcf"
  }

  runtime {
    docker: "getwilds/varscan:2.4.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}


task mpileup2cns {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run VarScan mpileup2cns for germline variant calling"
    outputs: {
        vcf: "VCF file containing SNV and indel calls"
    }
  }

  parameter_meta {
    sample_name: "Name of the sample (used in output file names)"
    pileup: "Samtools mpileup file generated with `--no-BAQ` and reference FASTA"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    String sample_name
    File pileup
    Int memory_gb = 16
    Int cpu_cores = 4
  }

  command <<<
    set -euo pipefail

    java -Xmx~{memory_gb - 2}g -jar /usr/local/bin/VarScan.jar \
      mpileup2cns \
      --variants \
      --output-vcf 1 \
      > "~{sample_name}.vcf" \
      < "~{pileup}"
  >>>

  output {
    File vcf = "~{sample_name}.vcf"
  }

  runtime {
    docker: "getwilds/varscan:2.4.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
