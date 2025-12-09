# ww-gdc Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading genomic data from the NCI Genomic Data Commons (GDC) using the GDC Data Transfer Tool. This module enables easy integration of TCGA and other cancer genomics data downloads into WDL workflows, supporting both controlled-access and open-access data.

## Overview

The GDC Data Transfer Tool provides efficient, resumable downloads of large genomic datasets from the NCI GDC. This module wraps the `gdc-client` command-line tool to enable seamless integration with WDL workflows.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: [ww-gdc.wdl](ww-gdc.wdl) - Contains task definitions for the module
- **Test workflow**: [testrun.wdl](testrun.wdl) - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `download_by_manifest`

Download files from GDC using a manifest file. This is the recommended approach for batch downloads.

**Inputs:**
- `manifest_file` (File, required): GDC manifest file containing file UUIDs and metadata (downloadable from GDC Data Portal)
- `token_file` (File, optional): Authentication token file for controlled-access data (downloadable from GDC user profile)
- `n_processes` (Int, default=8): Number of parallel download processes
- `retry_amount` (Int, default=5): Number of times to retry failed downloads
- `wait_time` (Int, default=5): Seconds to wait between retries
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `downloaded_files` (Array[File]): Array of downloaded data files from GDC
- `download_log` (File): Log file with download statistics and any errors

### `download_by_uuids`

Download files from GDC using file UUIDs. Useful when you have specific file IDs rather than a manifest.

**Inputs:**
- `file_uuids` (Array[String], required): Array of GDC file UUIDs to download
- `token_file` (File, optional): Authentication token file for controlled-access data (downloadable from GDC user profile)
- `n_processes` (Int, default=8): Number of parallel download processes
- `retry_amount` (Int, default=5): Number of times to retry failed downloads
- `wait_time` (Int, default=5): Seconds to wait between retries
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `downloaded_files` (Array[File]): Array of downloaded data files from GDC
- `download_log` (File): Log file with download statistics and any errors

## Usage as a Module

### Getting a GDC Authentication Token

For controlled-access data, you'll need an authentication token:

1. Go to the [GDC Data Portal](https://portal.gdc.cancer.gov/)
2. Log in with your eRA Commons credentials
3. Click your username in the top right
4. Select "Download Token" from the dropdown menu
5. Save the token file securely (it's valid for 30 days)

**Important**: Keep your token file secure. It provides access to all controlled-access data you're authorized to view.

### Obtaining a Manifest File

To download data using a manifest:

1. Go to the [GDC Data Portal](https://portal.gdc.cancer.gov/)
2. Use the Repository or Exploration page to find your files of interest
3. Add files to your cart
4. In the cart, click "Manifest" to download the manifest file
5. The manifest is a tab-separated file containing file UUIDs, filenames, and metadata

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gdc/ww-gdc.wdl" as gdc

workflow my_tcga_analysis {
  input {
    File gdc_manifest
    File gdc_token
  }

  # Download TCGA data from GDC
  call gdc.download_by_manifest {
    input:
      manifest_file = gdc_manifest,
      token_file = gdc_token,
      n_processes = 16
  }

  # Process downloaded files in your workflow
  scatter (data_file in download_by_manifest.downloaded_files) {
    call my_analysis_task {
      input: input_file = data_file
    }
  }

  output {
    Array[File] analysis_results = my_analysis_task.results
  }
}
```

### Advanced Usage Examples

**Download open-access data without a token:**
```wdl
call gdc.download_by_manifest {
  input:
    manifest_file = open_access_manifest
    # No token_file needed for open-access data
}
```

**Download specific files by UUID:**
```wdl
call gdc.download_by_uuids {
  input:
    file_uuids = [
      "e4da3d4a-9f3c-4a0e-9a4e-c8f8a9b4d5e6",
      "f5eb4e5b-0g4d-5b1f-0b5f-d9g9b0c5e6f7"
    ],
    token_file = gdc_token
}
```

**Custom resource allocation for large downloads:**
```wdl
call gdc.download_by_manifest {
  input:
    manifest_file = large_manifest,
    token_file = gdc_token,
    n_processes = 32,
    cpu_cores = 16,
    memory_gb = 32
}
```

## Implementation Details

### Cromwell Path Length Workaround

This module includes a workaround for a Cromwell-specific limitation where deeply nested execution directories can exceed the Unix socket path length limit (~108 characters). The gdc-client tool uses Python multiprocessing which creates Unix domain sockets, and Cromwell's deep directory structures can cause "AF_UNIX path too long" errors.

**How the workaround works:**

1. Downloads execute in a short temporary path (`/tmp/gdc_$$`) instead of the Cromwell execution directory
2. Environment variables (`TMPDIR`, `TEMP`, `TMP`) are set to ensure Python multiprocessing creates sockets in the short path
3. Downloaded files are moved back to the Cromwell execution directory after download completes
4. Files are organized under `manifest/` or `uuids/` subdirectories to preserve the UUID-based directory structure

**Impact on users:**
- No changes to task inputs or outputs required
- Works seamlessly across miniWDL, Cromwell, and Sprocket
- Downloaded files maintain the same directory structure as native gdc-client

## Testing the Module

The module includes a test workflow ([testrun.wdl](testrun.wdl)) that demonstrates downloading open-access data:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint gdc_client_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

**Note**: The test workflow downloads a small open-access file and does not require authentication.

## Docker Container

This module uses the WILDS Docker Library gdc-client image, which includes:
- GDC Data Transfer Tool (gdc-client)
- Python runtime and dependencies
- All necessary system libraries

## Data Types Supported

The GDC Data Transfer Tool can download any file type available in the GDC, including:

- **Sequencing Data**: BAM, FASTQ, VCF files
- **Clinical Data**: XML, TSV, JSON files
- **Biospecimen Data**: Metadata files
- **Analysis Results**: MAF (mutation annotation format), segmentation files
- **Images**: Slide images for pathology data

## Performance Considerations

### Default Resources
- **CPU**: 4 cores
- **Memory**: 8 GB
- **Parallel Downloads**: 8 simultaneous connections
- **Runtime**: Varies based on file sizes and network speed

### Resource Scaling Recommendations

For large downloads (>100GB total):
- Increase `n_processes` to 16-32 for faster parallel downloads
- Increase `cpu_cores` to 8-16 to support more parallel processing
- Increase `memory_gb` to 16-32 GB for handling large manifests

For slow or unstable network connections:
- Increase `retry_amount` to 10-20
- Increase `wait_time` to 10-30 seconds
- Decrease `n_processes` to 4-8 to reduce simultaneous connections

## Troubleshooting

### Common Issues

**"401 Unauthorized" errors:**
- Your token file may be expired (tokens are valid for 30 days)
- Download a new token from the GDC Data Portal
- Ensure you're authorized to access the requested controlled-access data

**Download failures or incomplete transfers:**
- The GDC Data Transfer Tool automatically resumes interrupted downloads
- Check the `download.log` file for specific error messages
- Increase `retry_amount` and `wait_time` parameters

**"No data files found" errors:**
- Verify your manifest file is properly formatted
- Check that the file UUIDs in your manifest are valid
- Ensure the files haven't been removed from GDC

**"AF_UNIX path too long" errors (legacy issue):**
- This issue has been resolved in the current version through the Cromwell path workaround
- See the "Implementation Details" section above for technical details

## Citation

### GDC Data Transfer Tool

If you use this module to download data from the GDC, please cite:

> NCI Genomic Data Commons (GDC)
> National Cancer Institute
> https://gdc.cancer.gov

### TCGA Data

If you use TCGA data downloaded via this tool, please cite:

> The Cancer Genome Atlas (TCGA)
> National Cancer Institute
> https://www.cancer.gov/tcga

Please also follow the [GDC Data Use Policies](https://gdc.cancer.gov/access-data/data-access-policies) and cite the specific TCGA publication(s) relevant to your cancer type(s).

## Additional Resources

- **[GDC Data Portal](https://portal.gdc.cancer.gov/)**: Web interface for browsing and selecting data
- **[GDC Documentation](https://docs.gdc.cancer.gov/)**: Complete documentation for the GDC
- **[GDC Data Transfer Tool User Guide](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/)**: Detailed command-line usage
- **[GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)**: Programmatic data access
- **[GDC Data Access Policies](https://gdc.cancer.gov/access-data/data-access-policies)**: Terms of use and citation requirements

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

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
