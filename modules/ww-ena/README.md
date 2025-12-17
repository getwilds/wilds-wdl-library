# ww-ena Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading sequencing data files from the European Nucleotide Archive (ENA) using the [ena-file-downloader tool](https://github.com/enasequence/ena-ftp-downloader). This module enables automated retrieval of FASTQ, BAM, and other sequencing data formats directly from ENA using accession numbers or search queries.

## Overview

The European Nucleotide Archive (ENA) is one of the world's largest repositories of nucleotide sequence data. This module wraps the ena-file-downloader tool to provide seamless integration of ENA data downloads into WDL workflows, supporting both FTP and Aspera transfer protocols.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-ena.wdl` - Contains task definitions for downloading files from ENA
- **Test workflow**: `testrun.wdl` - Demonstration workflow with validation
- **Container**: `getwilds/ena-tools:2.1.1` from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)

## Available Tasks

### `download_files`

Downloads sequencing data files from ENA using accession numbers. Supports single or multiple accessions provided as a comma-separated string or via a file.

**Inputs:**
- `accessions` (String, optional): Comma-separated list of ENA accession numbers (e.g., "ERR2208926,ERR2208890") or a single accession.
- `accessions_file` (File, optional): File containing accession numbers (one per line or tab-separated with accessions in first column)
- `file_format` (String, default="READS_FASTQ"): Format of files to download. Options are:
  - `READS_FASTQ`: FASTQ format reads
  - `READS_SUBMITTED`: Submitted read files
  - `READS_BAM`: BAM format alignments
  - `ANALYSIS_SUBMITTED`: Submitted analysis files
  - `ANALYSIS_GENERATED`: Generated analysis files
- `protocol` (String, default="FTP"): Transfer protocol (`FTP` or `ASPERA`)
- `aspera_location` (String, optional): Path to Aspera Connect/CLI installation (required if protocol is ASPERA)
- `output_dir_name` (String, default="ena_downloads"): Name for the output directory
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Note:** Either `accessions` or `accessions_file` must be provided.

**Outputs:**
- `downloaded_files` (Array[File]): Array of downloaded files from ENA
- `download_log` (File): Log file containing download status and details
- `download_summary` (File): Summary report of the download operation
- `accessions_used` (String): The accession numbers that were processed

### `download_by_query`

Downloads sequencing data files from ENA using a search query. Allows filtering by result type and query parameters to retrieve multiple files matching specific criteria.

**Inputs:**
- `query` (String): ENA search query string containing result type and query parameters (e.g., 'result=read_run&query=study_accession=PRJEB1234')
- `file_format` (String, default="READS_FASTQ"): Format of files to download (same options as `download_files`)
- `protocol` (String, default="FTP"): Transfer protocol (`FTP` or `ASPERA`)
- `aspera_location` (String, optional): Path to Aspera Connect/CLI installation
- `output_dir_name` (String, default="ena_downloads"): Name for the output directory
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `downloaded_files` (Array[File]): Array of downloaded files from ENA
- `download_log` (File): Log file containing download status and details
- `download_summary` (File): Summary report of the download operation

### `extract_fastq_pairs`

Extracts R1 and R2 FASTQ files from downloaded ENA files for downstream paired-end processing. This task identifies paired-end FASTQ files by common naming patterns, creates standardized outputs, and automatically extracts the accession ID from the filename.

**Inputs:**
- `downloaded_files` (Array[File]): Array of files downloaded from ENA (typically from `download_files` task)

**Outputs:**
- `r1` (File): Read 1 FASTQ file
- `r2` (File): Read 2 FASTQ file
- `accession` (String): ENA accession ID automatically extracted from the filename

**Usage Note:** This task is designed for FASTQ workflows requiring separate R1/R2 files. It searches for common paired-end naming patterns including `_1.fastq.gz`/`_2.fastq.gz`, `_R1.fastq.gz`/`_R2.fastq.gz`, and their uncompressed equivalents. The accession ID is automatically extracted from the filename (e.g., `ERR000001_1.fastq.gz` → `ERR000001`). If you're downloading other file formats (BAM, analysis files), you don't need this task.

## Usage as a Module

### Basic Usage: Download by Accession Numbers

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ena_tasks

workflow download_ena_data {
  input {
    Array[String] sample_accessions
  }

  # Download FASTQ files for each accession
  scatter (accession in sample_accessions) {
    call ena_tasks.download_files {
      input:
        accessions = accession,
        file_format = "READS_FASTQ",
        protocol = "FTP"
    }
  }

  output {
    Array[Array[File]] all_downloads = download_files.downloaded_files
  }
}
```

### Download Multiple Accessions from a File

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ena_tasks

workflow download_from_file {
  input {
    File accession_list
  }

  call ena_tasks.download_files {
    input:
      accessions_file = accession_list,
      file_format = "READS_FASTQ",
      output_dir_name = "my_ena_downloads"
  }

  output {
    Array[File] fastq_files = download_files.downloaded_files
    File summary = download_files.download_summary
  }
}
```

### Download Using a Search Query

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ena_tasks

workflow download_by_study {
  input {
    String study_accession
  }

  call ena_tasks.download_by_query {
    input:
      query = "result=read_run&query=study_accession=~{study_accession}",
      file_format = "READS_FASTQ"
  }

  output {
    Array[File] study_fastqs = download_by_query.downloaded_files
  }
}
```

### Extract FASTQ Pairs for Paired-End Analysis

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ena_tasks

workflow ena_paired_end_workflow {
  input {
    Array[String] sample_accessions
  }

  scatter (accession in sample_accessions) {
    # Download FASTQ files from ENA
    call ena_tasks.download_files {
      input:
        accessions = accession,
        file_format = "READS_FASTQ"
    }

    # Extract R1 and R2 for paired-end processing
    # Accession ID is automatically extracted from filenames
    call ena_tasks.extract_fastq_pairs {
      input:
        downloaded_files = download_files.downloaded_files
    }
  }

  output {
    Array[File] all_r1 = extract_fastq_pairs.r1
    Array[File] all_r2 = extract_fastq_pairs.r2
    Array[String] all_accessions = extract_fastq_pairs.accession
  }
}
```

### Integration with Downstream Analysis

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ena_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastqc/ww-fastqc.wdl" as fastqc_tasks

workflow ena_to_qc_pipeline {
  input {
    String ena_accession
  }

  # Download data from ENA
  call ena_tasks.download_files {
    input:
      accessions = ena_accession
  }

  # Run quality control on downloaded files
  scatter (fastq_file in download_files.downloaded_files) {
    call fastqc_tasks.fastqc {
      input:
        fastq = fastq_file
    }
  }

  output {
    Array[File] qc_reports = fastqc.html_report
  }
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint ena_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Validation

The test workflow includes a validation task that checks:
- Files were successfully downloaded
- Output files are non-empty
- Download logs and summaries are generated

## Docker Container

This module uses the `getwilds/ena-tools:2.1.1` container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- ena-file-downloader tool (version 2.1.1)
- Java runtime environment
- Bash shell
- Standard Unix utilities

## ENA Data Access

### Finding Accession Numbers

ENA accession numbers can be found through:
- **ENA Browser**: https://www.ebi.ac.uk/ena/browser/
- **ENA Search**: https://www.ebi.ac.uk/ena/browser/advanced-search

Common accession formats:
- Study: `PRJEB*`, `PRJNA*`, `ERP*`, `SRP*`
- Sample: `SAMEA*`, `SAMN*`, `ERS*`, `SRS*`
- Run: `ERR*`, `SRR*`, `DRR*`

### File Formats

The module supports multiple data formats:
- **READS_FASTQ**: Most common format for sequencing reads
- **READS_SUBMITTED**: Original submitted format (may be FASTQ, FASTA, or other)
- **READS_BAM**: Aligned reads in BAM format
- **ANALYSIS_SUBMITTED**: Analysis files as submitted
- **ANALYSIS_GENERATED**: ENA-generated analysis files

## Performance Considerations

### Transfer Protocols

- **FTP (default)**: Reliable, works everywhere, moderate speed
- **ASPERA**: Much faster for large files, requires Aspera Connect/CLI installation

### Resource Requirements

Default settings work well for most use cases:
- **CPU**: 2 cores (downloads are not CPU-intensive)
- **Memory**: 8 GB (sufficient for download management)
- **Disk**: Varies by data size; ensure adequate storage for your downloads

### Runtime Estimates

Download times depend on:
- File sizes (FASTQ files typically 1-10 GB per sample)
- Network speed
- Transfer protocol (Aspera is significantly faster)
- ENA server load

Typical download speeds:
- **FTP**: 10-50 MB/s
- **ASPERA**: 100-500 MB/s (with good network connection)

## Citation

If you use this module in your research, please cite:

**ENA Database:**
> Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration.
> The sequence read archive. Nucleic Acids Res. 2011 Jan;39(Database issue):D19-21.
> DOI: 10.1093/nar/gkq1019

**ena-file-downloader tool:**
> https://github.com/enasequence/ena-ftp-downloader

## Common Use Cases

### Download Paired-End FASTQ Files

```wdl
call ena_tasks.download_files {
  input:
    accessions = "ERR2208926",  # This will download both _1 and _2 files
    file_format = "READS_FASTQ"
}
```

### Download Entire Study

```wdl
call ena_tasks.download_by_query {
  input:
    query = "result=read_run&query=study_accession=PRJEB1234",
    file_format = "READS_FASTQ"
}
```

### Download with Accession List File

Create a file `accessions.txt`:
```
ERR2208926
ERR2208890
ERR2208960
```

```wdl
call ena_tasks.download_files {
  input:
    accessions_file = "accessions.txt",
    file_format = "READS_FASTQ"
}
```

## Troubleshooting

### Download Failures

If downloads fail, check:
1. Accession numbers are valid (verify on ENA website)
2. Network connectivity to ENA servers
3. Adequate disk space for downloads
4. Review `download.log` for specific error messages

### Missing Files

Some accessions may not have all data types available. Check:
- The ENA browser for available data formats
- Try `READS_SUBMITTED` if `READS_FASTQ` is unavailable
- Review the download summary for details

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[ENA Documentation](https://www.ebi.ac.uk/ena/browser/about)**: Official ENA documentation
- **[ena-file-downloader GitHub](https://github.com/enasequence/ena-ftp-downloader)**: Source code and documentation for the downloader tool
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
