## WILDS WDL module for variant annotation using Annovar.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task annovar_annotate {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Annotate variants using Annovar with customizable protocols and operations"
    outputs: {
        annotated_vcf: "VCF file with Annovar annotations added",
        annotated_table: "Tab-delimited table with variant annotations"
    }
  }

  parameter_meta {
    vcf_to_annotate: "Input VCF file to be annotated"
    ref_name: "Reference genome build name for Annovar annotation (e.g., 'hg38', 'hg19')"
    annovar_protocols: "Comma-separated list of annotation protocols to apply"
    annovar_operation: "Comma-separated list of operations corresponding to the protocols"
    cpu_cores: "Number of CPU cores to allocate for the annotation task"
    memory_gb: "Memory in GB to allocate for the annotation task"
  }

  input {
    File vcf_to_annotate
    String ref_name
    String annovar_protocols
    String annovar_operation
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  String base_vcf_name = basename(vcf_to_annotate, ".vcf.gz")

  command <<<
    set -eo pipefail
    perl /annovar/table_annovar.pl "~{vcf_to_annotate}" /annovar/humandb/ \
      -buildver "~{ref_name}" \
      -outfile "~{base_vcf_name}" \
      -remove \
      -protocol "~{annovar_protocols}" \
      -operation "~{annovar_operation}" \
      -nastring . -vcfinput
    sed -i "s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10\tOtherinfo11\tOtherinfo12\tOtherinfo13/Otherinfo/g" "~{base_vcf_name}.~{ref_name}_multianno.txt"
  >>>

  output {
    File annotated_vcf = "~{base_vcf_name}.~{ref_name}_multianno.vcf"
    File annotated_table = "~{base_vcf_name}.~{ref_name}_multianno.txt"
  }

  runtime {
    docker: "getwilds/annovar:~{ref_name}"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
