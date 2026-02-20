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
        "fastq_to_cram.batch_info": format_inputs_batch(ds),
        "fastq_to_cram.cpu_cores": ds.params["cpu_cores"],
        "fastq_to_cram.memory_gb": ds.params["memory_gb"]
    }
    write_json("inputs.0.json", inputs)


def format_inputs_batch(ds: PreprocessDataset):
    """
    Format the inputs using the SampleData FastqGroup structs. Note that 
    library_name and sequencing_center, if provided, will be the same for
    all FastqGroups. This is because those in puts are provided via the 
    web form and cannot be sample-specific.

    [
        {
            "sample_name": "SampleA",
            "library_name": "Lib_A",
            "sequencing_center": "SeqCenter",
            "fastq_groups": [
                {
                    "group_name": "SampleA",
                    "fastq_r1_locations": ["/path/to/SampleA_R1.fastq.gz"],
                    "fastq_r2_locations": ["/path/to/SampleA_R2.fastq.gz"]
                }
            ]
        },
        {
            "sample_name": "SampleB",
            "library_name": "Lib_A",
            "sequencing_center": "SeqCenter",
            "fastq_groups": [
                {
                    "group_name": "SampleB",
                    "fastq_r1_locations": ["/path/to/SampleB_R1.fastq.gz"],
                    "fastq_r2_locations": ["/path/to/SampleB_R2.fastq.gz"]
                }
            ]
        }
    ]
    """

    # Get a table listing the FASTQ files selected by the user.
    # This DataFrame has columns: sample, fastq_1, fastq_2
    df = ds.pivot_samplesheet()

    library_name = ds.params.get("library_name")
    sequencing_center = ds.params.get("sequencing_center")

    batch_info = []
    for _, r in df.iterrows():
        fastq_r1 = [r["fastq_1"]]
        fastq_r2 = [r["fastq_2"]]

        sample = {
            "sample_name": r["sample"],
            "fastq_groups": [
                {
                    "group_name": r["sample"],
                    "fastq_r1_locations": fastq_r1,
                    "fastq_r2_locations": fastq_r2
                }
            ]
        }
        if library_name:
            sample["library_name"] = library_name
        if sequencing_center:
            sample["sequencing_center"] = sequencing_center
        batch_info.append(sample)

    return batch_info


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
