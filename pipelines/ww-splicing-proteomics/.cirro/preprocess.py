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

    # Parse metadata to determine group assignments for each sample
    group_vals = parse_groups(ds)

    ds.logger.info("Formatting inputs")
    inputs = {
        "splicing_proteomics.samples": format_inputs_samples(ds, group_vals),
        "splicing_proteomics.reference_genome": format_inputs_reference_genome(ds),
        "splicing_proteomics.ensembl_gtf": ds.params["ensembl_gtf"],
        "splicing_proteomics.read_length": ds.params["read_length"],
        "splicing_proteomics.read_type": ds.params["read_type"],
        "splicing_proteomics.library_type": ds.params["library_type"],
        "splicing_proteomics.output_prefix": ds.params["output_prefix"],
        "splicing_proteomics.sjdb_overhang": ds.params["sjdb_overhang"],
        "splicing_proteomics.genome_sa_index_nbases": ds.params["genome_sa_index_nbases"],
        "splicing_proteomics.rmats_anchor_length": ds.params["rmats_anchor_length"],
        "splicing_proteomics.rmats_novel_splice_sites": ds.params["rmats_novel_splice_sites"],
        "splicing_proteomics.rmats_cstat": ds.params["rmats_cstat"],
        "splicing_proteomics.jcast_min_read_count": ds.params["jcast_min_read_count"],
        "splicing_proteomics.jcast_qvalue_min": ds.params["jcast_qvalue_min"],
        "splicing_proteomics.jcast_qvalue_max": ds.params["jcast_qvalue_max"],
        "splicing_proteomics.star_cpu": ds.params["star_cpu"],
        "splicing_proteomics.star_memory_gb": ds.params["star_memory_gb"],
        "splicing_proteomics.rmats_cpu": ds.params["rmats_cpu"],
        "splicing_proteomics.rmats_memory_gb": ds.params["rmats_memory_gb"],
        "splicing_proteomics.jcast_cpu": ds.params["jcast_cpu"],
        "splicing_proteomics.jcast_memory_gb": ds.params["jcast_memory_gb"],
    }

    # Only add jcast_splice_types if user specified specific types
    if ds.params.get("jcast_splice_types"):
        inputs["splicing_proteomics.jcast_splice_types"] = ds.params["jcast_splice_types"]

    write_json("inputs.0.json", inputs)


def parse_groups(ds: PreprocessDataset) -> dict:
    """
    The user provided metadata for each sample with a column indicating group
    assignment. Validate that the specified column exists and contains exactly
    two groups: "group1" and "group2".
    """

    # Get the metadata provided by the user as a pandas DataFrame
    meta: pd.DataFrame = ds.samplesheet
    ds.logger.info(f"User-provided metadata:\n{meta.to_csv()}")

    # Make sure that the group column specified by the user is present
    group_cname = ds.params["group_column"]
    if group_cname not in meta.columns.values:
        raise ValueError(f"Group column not found in metadata: {group_cname}")

    # Validate that the column contains valid group values
    valid_groups = {"group1", "group2"}
    actual_groups = set(meta[group_cname].unique())
    if not actual_groups.issubset(valid_groups):
        invalid = actual_groups - valid_groups
        raise ValueError(
            f"Invalid group values in '{group_cname}' column: {invalid}. "
            f"Values must be 'group1' or 'group2'."
        )
    if len(actual_groups) < 2:
        raise ValueError(
            f"Both 'group1' and 'group2' must be present in '{group_cname}' column. "
            f"Found only: {actual_groups}"
        )

    # Return a dictionary mapping sample name to group
    return meta.set_index("sample")[group_cname].to_dict()


def format_inputs_samples(ds: PreprocessDataset, group_vals: dict):
    """
    Format the inputs using the SampleInfo struct syntax:

    [
        {
            "name": "SRR1039508",
            "r1": "/path/to/SRR1039508_1.fastq.gz",
            "r2": "/path/to/SRR1039508_2.fastq.gz",
            "group": "group1"
        },
        ...
    ]
    """

    # Get a table listing the FASTQ files selected by the user
    # This DataFrame has columns: sample, fastq_1, fastq_2
    df = ds.pivot_samplesheet()

    return [
        {
            "name": r["sample"],
            "r1": r["fastq_1"],
            "r2": r["fastq_2"],
            "group": group_vals[r["sample"]],
        }
        for _, r in df.iterrows()
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
