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
        "imputation.chromosomes": format_chromosomes(ds),
        "imputation.input_crams": [v.strip() for v in ds.params["input_crams"].split(",")],
        "imputation.input_cram_indices": [v.strip() for v in ds.params["input_cram_indices"].split(",")],
        "imputation.reference_fasta": ds.params["reference_fasta"],
        "imputation.reference_fasta_index": ds.params["reference_fasta_index"],
        "imputation.output_prefix": ds.params["output_prefix"],
        "imputation.output_format": ds.params["output_format"],
        "imputation.chunk_cpu_cores": ds.params["chunk_cpu_cores"],
        "imputation.chunk_memory_gb": ds.params["chunk_memory_gb"],
        "imputation.phase_cpu_cores": ds.params["phase_cpu_cores"],
        "imputation.phase_memory_gb": ds.params["phase_memory_gb"],
        "imputation.ligate_cpu_cores": ds.params["ligate_cpu_cores"],
        "imputation.ligate_memory_gb": ds.params["ligate_memory_gb"],
        "imputation.concat_cpu_cores": ds.params["concat_cpu_cores"],
        "imputation.concat_memory_gb": ds.params["concat_memory_gb"],
    }
    write_json("inputs.0.json", inputs)


def format_chromosomes(ds: PreprocessDataset):
    """
    Format the chromosomes input using the ChromosomeData struct syntax:

    [
        {
            "chromosome": "chr22",
            "reference_vcf": "/path/to/1000GP_chr22.bcf",
            "reference_vcf_index": "/path/to/1000GP_chr22.bcf.csi",
            "genetic_map": "/path/to/genetic_map_chr22.txt"
        },
        ...
    ]
    """

    # Use the ds.files internal samplesheet, since file order isn't maintained
    # from the user uplaoded sample sheet and we have to identify files by extension.
    df = ds.files

    chromosomes = []
    for sample, grp in df.groupby("sample"):
        files = grp["file"].tolist()
        chromosomes.append({
            "chromosome": sample,
            "reference_vcf": next(f for f in files if f.endswith(".bcf") and not f.endswith(".bcf.csi")),
            "reference_vcf_index": next(f for f in files if f.endswith(".bcf.csi")),
            "genetic_map": next(f for f in files if f.endswith(".gmap.gz")),
        })

    return chromosomes


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
