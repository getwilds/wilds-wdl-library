#!/usr/bin/env python3
import json
from typing import Tuple
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

    # Parse the metadata provided by the user and make sure that the values
    # provided in the form are valid.
    meta_vals, cname, ref, comp = parse_contrasts(ds)

    ds.logger.info("Formatting inputs")
    inputs = {
        "star_deseq2.samples": format_inputs_samples(ds, meta_vals, cname, ref, comp),
        "star_deseq2.reference_genome": format_inputs_reference_genome(ds),
        "star_deseq2.reference_level": ref,
        "star_deseq2.contrast": f"condition,{comp},{ref}",
        "star_deseq2.star_cpu": ds.params["star_cpu"],
        "star_deseq2.star_memory_gb": ds.params["star_memory_gb"]
    }
    write_json("inputs.0.json", inputs)


def parse_contrasts(ds: PreprocessDataset) -> Tuple[pd.DataFrame, str, str, str]:
    """
    The user provided metadata for each sample, and also filled out the form
    when submitting the form.
    Let's make sure that all of those selections are valid, and that they
    make sense together.
    """

    # Get the metadata provided by the user annotating these samples
    # as a pandas DataFrame.
    meta: pd.DataFrame = ds.samplesheet
    ds.logger.info(f"User-provided metadata:\n{meta.to_csv()}")

    # Make sure that the column specified by the user for contasting
    # samples is present.
    cname = ds.params["contrast_cname"]
    # Stop execution if it is not present and inform the user of the error
    if cname not in meta.columns.values:
        raise ValueError(f"Metadata Contrast Reference Level not found: {cname}")

    # Get the group names provided by the user for contrasting samples
    ref = ds.params["contrast_reference_level"]
    comp = ds.params["contrast_comparison_level"]
    assert ref != comp, "Reference and comparison groups cannot be the same"
    assert ref in meta[cname].values, f"Value not found in {cname} column: {ref}"
    assert comp in meta[cname].values, f"Value not found in {cname} column: {comp}"

    # Get the metadata annotations for each sample as a dictionary
    meta_vals = meta.set_index("sample")[cname].to_dict()

    return meta_vals, cname, ref, comp


def format_inputs_samples(
    ds: PreprocessDataset,
    meta_vals: pd.DataFrame,
    cname: str,
    ref: str,
    comp: str
):
    """
    Format the inputs using the syntax:

    [
        {
            "name": "SRR10724344",
            "r1": "/path_to_first_input/SRR10724344_1.fastq.gz",
            "r2": "/path_to_first_input/SRR10724344_2.fastq.gz",
            "condition": "treatment"
        },
        {
            "name": "SRR12345678",
            "r1": "/path_to_second_input/SRR12345678_1.fastq.gz",
            "r2": "/path_to_second_input/SRR12345678_2.fastq.gz",
            "condition": "control"
        }
    ]
    """

    # Get a table listing the FASTQ files selected by the user.
    # This DataFrame has columns: sample, fastq_1, fastq_2
    df = ds.pivot_samplesheet()

    # Reformat this as a list, only including the samples which belong
    # to the reference and comparison groups
    return [
        {
            "name": r["sample"],
            "r1": r["fastq_1"],
            "r2": r["fastq_2"],
            "condition": meta_vals[r["sample"]],
        }
        for _, r in df.iterrows()
        if meta_vals[r["sample"]] in [ref, comp]
    ]


def format_inputs_reference_genome(ds: PreprocessDataset):
    """
    Using the selection of the user, fill in the reference genome information
    """
    return {
        "name": "ref_genome",
        "fasta": ds.params["fasta"],
        "gtf": ds.params["gtf"]
    }


def setup_options(ds: PreprocessDataset):
    options = {
        "workflow_failure_mode": ds.params["workflow_failure_mode"],
        "write_to_cache": True,
        "read_from_cache": True,
        "default_runtime_attributes": {
            "maxRetries": 1
        },
        "final_workflow_outputs_dir": ds.params["final_workflow_outputs_dir"],
        "use_relative_output_paths": True
    }
    write_json("options.json", options)


if __name__ == "__main__":
    main()
