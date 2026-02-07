#!/usr/bin/env python3
import json
import logging
import os
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset, read_json


def main():
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)-8s [PreprocessDataset] %(message)s'
    )
    logger = logging.getLogger(__name__)

    # For workflows without input datasets, we can't use from_running()
    # because it expects a files.csv. Instead, we read params directly.
    dataset_root = os.getenv("PW_S3_DATASET")
    config_dir = f"{dataset_root}/config"

    logger.info(f"Reading params from {config_dir}")

    # Read params from the config directory using the SDK's S3-aware read_json
    params = read_json(f"{config_dir}/params.json")

    # Read metadata if it exists (wrap in try/except for S3)
    try:
        metadata = read_json(f"{config_dir}/metadata.json")
    except Exception:
        metadata = {}

    # Create PreprocessDataset with empty files/samplesheet since this
    # workflow downloads data from SRA rather than using input datasets
    ds = PreprocessDataset(
        samplesheet=pd.DataFrame(columns=["sample"]),
        files=pd.DataFrame(columns=["sample", "file"]),
        params=params,
        metadata=metadata,
        dataset_root=dataset_root
    )

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

    # Parse the comma-separated SRA IDs into a list
    sra_ids = [id.strip() for id in ds.params["sra_ids"].split(",")]
    ds.logger.info(f"Processing {len(sra_ids)} SRA samples: {sra_ids}")

    inputs = {
        "sra_star.sra_id_list": sra_ids,
        "sra_star.ref_genome": format_inputs_reference_genome(ds),
        "sra_star.sjdb_overhang": ds.params["sjdb_overhang"],
        "sra_star.genome_sa_index_nbases": ds.params["genome_sa_index_nbases"],
        "sra_star.ncpu": ds.params["ncpu"],
        "sra_star.memory_gb": ds.params["memory_gb"]
    }

    # Add max_reads only if it was specified (not None or empty)
    if ds.params.get("max_reads"):
        inputs["sra_star.max_reads"] = ds.params["max_reads"]

    write_json("inputs.0.json", inputs)


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
