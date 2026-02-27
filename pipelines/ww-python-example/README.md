# ww-python-example Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A minimal WILDS WDL pipeline demonstrating how to wrap a Python script in a WDL workflow, passing string, integer, and file inputs as command-line arguments.

## Overview

This pipeline provides a simple, self-contained example of executing a user-supplied Python script within a WDL task. It is intended as a starting point for collaborators who want to understand the basics of WDL and how to integrate their own scripts into a workflow.

The pipeline accepts a Python script along with three types of inputs (a string, an integer, and a file), passes them as command-line arguments to the script, and captures the output file.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

## Pipeline Steps

1. **Run Python Script** (`run_script` task):
   - Executes the provided Python script using `python3`
   - Passes the string input via `--name`
   - Passes the integer input via `--iterations`
   - Passes the file input via `--input-file`
   - Captures the output written to `results.txt`

## Key WDL Concepts Demonstrated

- **`File` inputs**: WDL localizes (copies) input files into the task's working directory, so your script receives a valid local path
- **`String` and `Int` inputs**: Primitive types are substituted directly into the command block
- **`~{variable}` syntax**: How WDL injects input values into the shell command within a `command <<<` block
- **`runtime` block**: Specifies the Docker image, CPU, and memory for the task
- **`output` block**: Declares which files the task produces for downstream use

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support

### Input Configuration

An example inputs file is provided in `inputs.json`. Create your own inputs JSON file with paths to your Python script and input data:

```json
{
  "run_python_script.python_script": "/path/to/example_script.py",
  "run_python_script.sample_name": "sample1",
  "run_python_script.num_iterations": 10,
  "run_python_script.input_file": "/path/to/input_data.txt"
}
```

An example Python script (`example_script.py`) is included to demonstrate the expected argument interface.

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-python-example.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-python-example.wdl -i inputs.json

# Using Sprocket
sprocket run ww-python-example.wdl inputs.json
```

### Running the Test Workflow

To test the pipeline without providing any external inputs:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow automatically generates a small Python script and input file inline, then runs the main workflow with them.

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `python_script` | The Python script to execute | File | Yes | - |
| `sample_name` | A string argument passed to the script via `--name` | String | Yes | - |
| `num_iterations` | An integer argument passed to the script via `--iterations` | Int | Yes | - |
| `input_file` | A file argument passed to the script via `--input-file` | File | Yes | - |

## Output Files

| Output | Description |
|--------|-------------|
| `results` | The output file (`results.txt`) produced by the Python script |

## Adapting This Example

To wrap your own Python script:

1. **Update your script** to accept `argparse`-style command-line arguments (see `example_script.py` for a template)
2. **Edit the `command <<<` block** in `ww-python-example.wdl` to match your script's argument names
3. **Add or remove inputs** in both the `workflow` and `task` `input` blocks as needed
4. **Update the `output` block** to match the filename(s) your script produces
5. **Change the Docker image** in the `runtime` block if your script requires additional Python packages

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
