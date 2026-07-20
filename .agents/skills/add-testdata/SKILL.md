---
name: add-testdata
description: Add a new test data download task to the ww-testdata module for use in module/pipeline testruns
argument-hint: <task-name>
---

# Add a New Test Data Task to ww-testdata

Add a new task called `$ARGUMENTS` to `modules/ww-testdata/ww-testdata.wdl`.

## Steps

### 1. Understand the Request

Clarify with the user if not obvious:
- What data needs to be downloaded or generated?
- Where is the source (public URL, S3 bucket, GitHub repo, generated synthetically)?
- What output files should the task produce?
- Which module or pipeline will use this test data?

### 2. Study Existing Patterns

Read `modules/ww-testdata/ww-testdata.wdl` to understand the existing task patterns. Key conventions:

**Docker images used for test data tasks:**
- `getwilds/samtools:1.11` — reference genome prep, BAM/CRAM handling, general file downloads via curl/wget
- `getwilds/awscli:2.27.49` — AWS S3 downloads (no-sign-request), general downloads
- `getwilds/bcftools:1.19` — VCF filtering, region extraction, normalization
- `getwilds/deseq2:1.40.2` — R/Bioconductor data generation
- `getwilds/r-utils:0.1.0` — general R scripts

**Common data sources:**
- UCSC (hgdownload.soe.ucsc.edu) — reference genomes
- GATK test data bucket (s3://gatk-test-data) — FASTQ, BAM
- Public GitHub repos — tool-specific test data
- 1000 Genomes FTP — population VCFs
- NCBI FTP — dbSNP
- UniProt — protein sequences
- Google Cloud Storage — GATK best practices

**Optimization strategies (keep test data small!):**
- Subsample BAM files (5-10% with `samtools view -s`)
- Filter to specific regions (`samtools faidx`, `bcftools view -r`)
- Use small chromosomes or chromosome subsets
- Download only what's needed (specific columns, specific samples)

### 3. Write the New Task

Add the task to `modules/ww-testdata/ww-testdata.wdl` following this template:

```wdl
task $ARGUMENTS {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Description of what this test data task provides"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
      output_name: "Description of each output file"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    # Add task-specific inputs with defaults where possible
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Download or generate test data
    # Keep files small for CI testing!
  >>>

  output {
    # Define typed outputs
  }

  runtime {
    docker: "getwilds/<appropriate-image>:<version>"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
```

### 4. Update the testrun.wdl

Read `modules/ww-testdata/testrun.wdl` and add a call to the new task. The testdata testrun calls ALL tasks to validate they work. Add:
- A call to the new task
- Include its outputs in the `validate_outputs` task if one exists

### 5. Update the README

Read `modules/ww-testdata/README.md` and add documentation for the new task, including:
- Task name and description
- Inputs (with defaults)
- Outputs (with types and descriptions)
- Data source

### 6. Lint

```bash
make lint NAME=ww-testdata
```

Fix any issues and re-run until clean.

### 7. Summary

Report:
- Task name and what it provides
- Data source and approximate size
- Which Docker image was used
- Which module(s)/pipeline(s) can now use this test data
