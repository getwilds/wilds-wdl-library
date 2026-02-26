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
        "ww_leukemia.samples": format_inputs_samples(ds),
        "ww_leukemia.ref_fasta": ds.params["ref_fasta"],
        "ww_leukemia.ref_fasta_index": ds.params["ref_fasta_index"],
        "ww_leukemia.ref_dict": ds.params["ref_dict"],
        "ww_leukemia.ref_name": ds.params["ref_name"],
        "ww_leukemia.dbsnp_vcf": ds.params["dbsnp_vcf"],
        "ww_leukemia.af_only_gnomad": ds.params["af_only_gnomad"],
        "ww_leukemia.known_indels_sites_vcfs": format_list(ds.params["known_indels_sites_vcfs"]),
        "ww_leukemia.wig_gc": ds.params["wig_gc"],
        "ww_leukemia.wig_map": ds.params["wig_map"],
        "ww_leukemia.panel_of_norm_rds": ds.params["panel_of_norm_rds"],
        "ww_leukemia.centromeres": ds.params["centromeres"],
        "ww_leukemia.ichorcna_chromosomes": ds.params["ichorcna_chromosomes"],
        "ww_leukemia.ichorcna_chrs_string": ds.params["ichorcna_chrs_string"],
        "ww_leukemia.annovar_protocols": ds.params["annovar_protocols"],
        "ww_leukemia.annovar_operation": ds.params["annovar_operation"],
        "ww_leukemia.scatter_count": ds.params["scatter_count"],
        "ww_leukemia.high_intensity_cpus": ds.params["high_intensity_cpus"],
        "ww_leukemia.high_intensity_memory_gb": ds.params["high_intensity_memory_gb"],
        "ww_leukemia.standard_cpus": ds.params["standard_cpus"],
        "ww_leukemia.standard_memory_gb": ds.params["standard_memory_gb"],
    }
    write_json("inputs.0.json", inputs)


def format_list(value: str) -> list:
    """Parse a comma-separated string into a list."""
    return [v.strip() for v in value.split(",")]


def format_inputs_samples(ds: PreprocessDataset):
    """
    Format the inputs using the SampleDetails struct syntax:

    [
        {
            "name": "SampleA",
            "cramfiles": ["/path/to/SampleA_run1.cram", "/path/to/SampleA_run2.cram"]
        },
        {
            "name": "SampleB",
            "cramfiles": ["/path/to/SampleB.cram"]
        }
    ]
    """

    # Get the samplesheet with CRAM files.
    # Expected columns: sample, file
    df = ds.samplesheet

    samples = []
    for sample_name, group in df.groupby("sample"):
        cram_files = group["file"].tolist()
        samples.append({
            "name": sample_name,
            "cramfiles": cram_files
        })

    return samples


def setup_options(ds: PreprocessDataset):
    options = {
        "workflow_failure_mode": ds.params["workflow_failure_mode"],
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
