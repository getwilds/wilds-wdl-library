# ww-aws-sso
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for AWS operations using traditional AWS credentials and SSO authentication via the [AWS CLI](https://aws.amazon.com/cli/).

## Overview

This module provides reusable WDL tasks for common AWS operations including downloading files from S3 buckets (both public and private), uploading files to S3, listing bucket contents, and syncing directories. The module supports both authenticated operations using traditional AWS credentials/SSO and public bucket access without authentication.

Designed to be a foundational component within the WILDS ecosystem, this module is suitable for data management workflows, preprocessing pipelines, and integration into larger bioinformatics analyses that require cloud storage operations.

Certificate-based authentication via IAM Roles Anywhere will be available soon via a companion module `ww-aws-iamra`.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `s3_download_file`, `s3_upload_file`, `s3_list_bucket`, `validate_outputs`
- **Workflow**: `aws_sso_example` (demonstration workflow executing all tasks)
- **Container**: `getwilds/awscli:2.27.49`

## Tasks

### `s3_download_file`
Downloads a single file from an S3 bucket with support for both public and private buckets.

**Inputs:**
- `s3_uri` (String): S3 URI of the file to download (e.g., s3://bucket/path/file.txt)
- `output_filename` (String?): Name for the downloaded file (optional, uses original name if not specified)
- `aws_config_file` (File?): Path to AWS config file (optional, uses --no-sign-request if not provided)
- `aws_credentials_file` (File?): Path to AWS credentials file (optional)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)
- `memory_gb` (Int): Memory allocation in GB (default: 2)

**Outputs:**
- `downloaded_file` (File): Downloaded file from S3

### `s3_upload_file`
Uploads a single file to an S3 bucket (requires AWS credentials).

**Inputs:**
- `file_to_upload` (File): File to upload to S3
- `s3_bucket` (String): S3 bucket name (without s3:// prefix)
- `s3_key` (String?): S3 key/path for the uploaded file (optional, uses filename if not specified)
- `aws_config_file` (File): Path to AWS config file (required for uploads)
- `aws_credentials_file` (File?): Path to AWS credentials file (optional)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)
- `memory_gb` (Int): Memory allocation in GB (default: 2)

**Outputs:**
- `s3_uri` (String): S3 URI of the uploaded file

### `s3_list_bucket`
Lists contents of an S3 bucket or prefix with comprehensive configuration options.

**Inputs:**
- `s3_uri` (String): S3 URI to list (e.g., s3://bucket/ or s3://bucket/prefix/)
- `aws_config_file` (File?): Path to AWS config file (optional, uses --no-sign-request if not provided)
- `aws_credentials_file` (File?): Path to AWS credentials file (optional)
- `recursive` (Boolean): List recursively (default: true)
- `human_readable` (Boolean): Use human-readable file sizes (default: true)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)
- `memory_gb` (Int): Memory allocation in GB (default: 2)

**Outputs:**
- `file_list` (File): Text file containing list of S3 objects
- `object_count` (Int): Number of objects found

### `validate_outputs`
Validates AWS operation outputs and generates a comprehensive summary report.

**Inputs:**
- `downloaded_files` (Array[File]?): Files that were downloaded from S3 (optional)
- `bucket_listing` (File?): File containing S3 bucket listing (optional)
- `object_count` (Int?): Number of objects found in bucket listing (optional)

**Outputs:**
- `report` (File): Validation report with AWS operation statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl" as aws_sso_tasks

struct S3File {
    String s3_uri
    String local_name
}

workflow my_data_processing_pipeline {
  input {
    Array[S3File] input_files
    String output_bucket
    File? aws_config_file
    File? aws_credentials_file
  }
  
  # Download files from S3
  scatter (file in input_files) {
    call aws_sso_tasks.s3_download_file {
      input:
        s3_uri = file.s3_uri,
        output_filename = file.local_name,
        aws_config_file = aws_config_file,
        aws_credentials_file = aws_credentials_file
    }
  }
  
  # Process files here (your custom analysis tasks)
  # ...
  
  # Upload results back to S3
  scatter (processed_file in processed_files) {
    call aws_sso_tasks.s3_upload_file {
      input:
        file_to_upload = processed_file,
        s3_bucket = output_bucket,
        aws_config_file = aws_config_file,
        aws_credentials_file = aws_credentials_file
    }
  }
  
  output {
    Array[File] downloaded_data = s3_download_file.downloaded_file
    Array[String] uploaded_uris = s3_upload_file.s3_uri
  }
}
```

### Advanced Usage Examples

**Public bucket access (no credentials required):**
```wdl
call aws_sso_tasks.s3_download_file {
  input:
    s3_uri = "s3://gatk-test-data/wgs_fastq/NA12878_20k/sample.fastq",
    # No aws_config_file provided - uses --no-sign-request automatically
}
```

**Private bucket with SSO credentials:**
```wdl
call aws_sso_tasks.s3_download_file {
  input:
    s3_uri = "s3://my-private-bucket/data/sample.fastq",
    aws_config_file = "/path/to/aws_config",
    aws_credentials_file = "/path/to/aws_credentials"
}
```

**Bucket inventory and selective downloading:**
```wdl
# First, list bucket contents
call aws_sso_tasks.s3_list_bucket {
  input:
    s3_uri = "s3://my-data-bucket/experiment-2024/",
    recursive = true,
    human_readable = true
}

# Then download specific files based on the listing
# (additional processing logic needed to parse file_list)
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **Data preprocessing**: Download reference genomes, test datasets, or input files
- **ww-sra**: Combine with SRA downloads for comprehensive data acquisition
- **ww-testdata**: Use alongside test data provisioning for development workflows
- **Result archival**: Upload analysis results, logs, and reports to S3 for long-term storage
- **Multi-stage pipelines**: Transfer intermediate results between workflow stages via S3

## Requirements

### **AWS Configuration**
- For public buckets: No configuration required (uses `--no-sign-request`)
- For private buckets: AWS credentials and configuration files must be mounted/available
- Sufficient permissions for the intended operations (read for downloads, write for uploads)

### **Infrastructure**
- Internet access for S3 operations
- Sufficient local storage for downloaded files
- Container runtime supporting the `getwilds/awscli` image

## Authentication Modes

The module automatically handles two authentication scenarios:

1. **Public Access**: When no `aws_config_file` is provided, operations use `--no-sign-request` for public bucket access
2. **Authenticated Access**: When `aws_config_file` is provided, standard AWS credential chain is used

### Example AWS Configuration Files

**AWS Config File (`~/.aws/config`):**
```ini
[default]
region = us-west-2
output = json
```

**AWS Credentials File (`~/.aws/credentials`):**
```ini
[default]
aws_access_key_id = AKIAIOSFODNN7EXAMPLE
aws_secret_access_key = wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
aws_session_token = IQoJb3JpZ2luX2VjEHAaCXVzLXdlc3QtMSJIMEYCIQDExample+TokenHereExampleToken==
```

**Note**: The `aws_session_token` is required when using temporary credentials (such as those obtained through AWS SSO, STS assume-role operations, or IAM roles). For permanent IAM user credentials, this field can be omitted.

## Security Considerations

### **Credential File Management**

**Important Security Notice**: When using AWS credentials with this module, copies of your credential and config files may be temporarily stored in the WDL execution engine's scratch directory. While these directories typically have restricted permissions, it's recommended to clean up your workspace after workflow completion to ensure credential files are completely removed.

**Best practices:**
- Delete workflow execution directories after completing your analysis
- Use IAM roles with minimal required permissions when possible
- Prefer temporary credentials over permanent access keys for enhanced security
- Regularly rotate AWS access keys used in workflows

**Not recommended:**
- Hardcoding credentials in WDL files
- Using root account credentials
- Sharing credentials between users or systems

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-aws-sso.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-aws-sso.wdl -i inputs.json

# Using Sprocket
sprocket run ww-aws-sso.wdl inputs.json
```

### Automatic Demo Mode

The demonstration workflow:
1. Lists contents of the GATK test data bucket (public access)
2. Downloads a sample FASTQ file from the public bucket
3. Validates all operations and generates a comprehensive report

## Performance Considerations

- **Network bandwidth**: Transfer speed depends on network connectivity and S3 region proximity
- **Concurrent operations**: Multiple parallel downloads/uploads scale with available CPU cores
- **Storage requirements**: Ensure adequate local disk space for downloaded files
- **Memory usage**: Scales conservatively with CPU allocation

## Error Handling

The module includes comprehensive error handling:
- Verification of successful file downloads
- Validation of upload operations
- Clear error messages for authentication issues
- Automatic detection of public vs. private bucket requirements

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real S3 operations with public test data
- Comprehensive validation of all outputs and error conditions
- Cross-platform testing to ensure container compatibility

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to AWS CLI usage, please refer to the [AWS CLI documentation](https://docs.aws.amazon.com/cli/).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
