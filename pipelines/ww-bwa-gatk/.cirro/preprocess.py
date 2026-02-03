#!/usr/bin/env python3
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def main():
    # Get the information provided by the user
    ds = PreprocessDataset.from_running()

    # Set up the inputs JSON
    setup_inputs(ds)
    # Set up the options JSON
    setup_options(ds)


def write_json(fp, obj, indent=4) -> None:

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def setup_inputs(ds: PreprocessDataset):
    """
    Format the inputs as JSON with the name "inputs.0.json"

    Note: The workflow could be invoked multiple times by writing out multiple
          inputs.*.json files. But for this scenario, only a single execution
          is needed.
    """

    ds.logger.info("Formatting inputs")

    inputs = {
        "bwa_gatk.samples": format_inputs_samples(ds),
        "bwa_gatk.reference_fasta": ds.params["reference_fasta"],
        "bwa_gatk.reference_fasta_index": ds.params["reference_fasta_index"],
        "bwa_gatk.dbsnp_vcf": ds.params["dbsnp_vcf"],
        "bwa_gatk.known_indels_vcf": ds.params["known_indels_vcf"],
        "bwa_gatk.cpu_cores": ds.params["cpu_cores"],
        "bwa_gatk.memory_gb": ds.params["memory_gb"]
    }
    write_json("inputs.0.json", inputs)


def format_inputs_samples(ds: PreprocessDataset):
    """
    Format the inputs using the BwaSample struct syntax:

    [
        {
            "name": "SRR10724344",
            "reads": "/path_to_first_input/SRR10724344_1.fastq.gz",
            "mates": "/path_to_first_input/SRR10724344_2.fastq.gz"
        },
        {
            "name": "SRR12345678",
            "reads": "/path_to_second_input/SRR12345678_1.fastq.gz"
        }
    ]
    """

    # Get a table listing the FASTQ files selected by the user.
    # This DataFrame has columns: sample, fastq_1, fastq_2
    df = ds.pivot_samplesheet()

    samples = []
    for _, r in df.iterrows():
        sample = {
            "name": r["sample"],
            "reads": r["fastq_1"],
        }
        # Only include mates if fastq_2 is present (paired-end)
        if "fastq_2" in r and pd.notna(r["fastq_2"]):
            sample["mates"] = r["fastq_2"]
        samples.append(sample)

    return samples


def setup_options(ds: PreprocessDataset):
    options = {
        "workflow_failure_mode": "ContinueWhilePossible",
        "write_to_cache": True,
        "read_from_cache": True,
        "default_runtime_attributes": {
            "maxRetries": 1
        },
        "final_workflow_outputs_dir": ds.params["out_dir"],
        "use_relative_output_paths": True
    }
    write_json("options.json", options)


if __name__ == "__main__":
    main()
