---
name: create-pipeline
description: Create a new WILDS WDL pipeline that combines existing modules into a workflow
argument-hint: <pipeline-name> <module1> <module2> [...]
allowed-tools: Bash, Edit, Read, Write, Glob, Grep, WebSearch, WebFetch
---

# Create a New WILDS WDL Pipeline

Create a new pipeline called `ww-$0` in `pipelines/ww-$0/` that combines the specified modules.

The remaining arguments (`$1`, `$2`, etc.) are the module names to combine. If no modules are specified, ask the user which modules to include.

## Steps

### 1. Verify Modules Exist

Check that each specified module exists under `modules/`. List available modules if any are missing:
```
modules/ww-<name>/ww-<name>.wdl
```

Read each module's WDL file to understand its tasks, inputs, and outputs. This is critical for wiring them together correctly.

### 2. Study Existing Pipelines

Read example pipelines to match the exact conventions. Study at least two of these:
- `pipelines/ww-sra-star/` (Basic: 2 modules)
- `pipelines/ww-star-deseq2/` (Intermediate: 4 modules)
- `pipelines/ww-bwa-gatk/` (Advanced: multi-step DNA-seq)

Pick the example closest in complexity to what you're building.

### 3. Design the Pipeline

Before writing code, plan:
- **Data flow:** How outputs from one module feed as inputs to the next
- **Scatter strategy:** What to scatter over (samples? chromosomes?) and what runs once (index building, reference prep)
- **Structs:** Group related inputs (e.g., `SampleInfo`, `RefGenome`) — follow existing struct patterns
- **Complexity level:** Basic (2-3 modules), Intermediate (4-6), or Advanced (10+)

Present the design to the user for confirmation before proceeding.

### 4. Create Pipeline Files

Create four files in `pipelines/ww-$0/`:

#### `ww-$0.wdl`
- `version 1.0`
- Import each module using GitHub raw URLs:
  ```wdl
  import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-<name>/ww-<name>.wdl" as <alias>_tasks
  ```
- Define structs for grouped inputs (samples, references)
- Workflow name: use the pipeline name with hyphens converted to underscores (e.g., `bwa_gatk`)
- Workflow MUST have a `meta` block with: author, email, description, url, outputs (same shape as task `meta` blocks in modules)
  - Use `author: "WILDS Team"` and `email: "wilds@fredhutch.org"` as placeholders — the actual contributor will replace these in review
- Wire module task outputs to downstream task inputs
- Include `cpu_cores` and `memory_gb` as workflow-level inputs with sensible defaults
- Collect final outputs at workflow level

#### `testrun.wdl`
- `version 1.0`
- Import the pipeline WDL and `ww-testdata`
- Workflow name: `{pipeline_name}_example` (underscored)
- Must run with zero configuration
- Use `ww-testdata` tasks for all test data — check what's available:
  - `download_ref_data` — reference genome, GTF, BED, index, dict
  - `download_fastq_data` — paired-end FASTQ
  - `download_bam_data` / `download_cram_data` — alignment files
  - `download_dbsnp_vcf` / `download_known_indels_vcf` / `download_gnomad_vcf` — variant databases
  - `download_test_transcriptome` — RNA transcriptome
  - `generate_pasilla_counts` — DESeq2 count data
  - `download_diamond_data` — protein sequence data
  - Various GLIMPSE2, ShapeMapper, CellRanger, ichorCNA tasks
- Use reduced resources for testing (fewer CPUs, less memory)
- Subsample data where possible (`max_reads`, region filtering)
- Construct structs from downloaded test data

#### `inputs.json`
- Use workflow name as prefix: `"pipeline_name.input_name"`
- Include placeholder paths: `"/path/to/file.ext"`
- Show struct format for complex inputs (arrays of structs)
- Include all configurable parameters with reasonable defaults
- Serve as a user template — make it clear what needs to be customized

#### `README.md`
- Follow existing pipeline README format (see `pipelines/ww-sra-star/README.md` or `pipelines/ww-star-deseq2/README.md`)
- Include: badges, overview, complexity level, pipeline structure, module dependencies, usage, input parameters, output files, testing section, citations

### 5. Lint the Pipeline

```bash
make lint NAME=ww-$0
```

Fix any issues and re-run until clean.

### 6. Test Run (if possible)

```bash
make run_sprocket NAME=ww-$0
```

### 7. Summary

Report:
- Pipeline complexity level and modules used
- Data flow diagram (ASCII)
- Files created
- Any issues needing follow-up (missing test data, Docker images, etc.)
