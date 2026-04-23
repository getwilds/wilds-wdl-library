#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :   compile_deseq2_results.py
@Time    :   2026/04/22
@Author  :   Taylor Firman
@Version :   0.1.0
@Contact :   tfirman@fredhutch.org
@Desc    :   Merge DESeq2 results with normalized counts and GTF gene annotations
"""

import pandas as pd
import re
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge DESeq2 results, normalized counts, and GTF annotations"
    )
    parser.add_argument(
        "--results", required=True, help="DESeq2 all-genes results CSV"
    )
    parser.add_argument(
        "--counts", required=True, help="DESeq2 normalized counts CSV"
    )
    parser.add_argument(
        "--gtf", required=True, help="GTF annotation file"
    )
    parser.add_argument(
        "--output",
        default="deseq2_compiled_results.csv",
        help="Output compiled CSV file",
    )
    return parser.parse_args()


def extract_attr(attributes, key):
    """Extract a specific attribute value from a GTF attributes string."""
    pattern = key + r'\s+"?([^";]+)"?'
    match = re.search(pattern, attributes)
    if match:
        return match.group(1)
    return None


def parse_gtf_annotations(gtf_path):
    """Parse gene annotations from a GTF file."""
    records = []
    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            record = {
                "gene_id": extract_attr(attrs, "gene_id"),
                "gene_name": extract_attr(attrs, "gene_name"),
                "gene_biotype": extract_attr(attrs, "gene_biotype"),
                "product": extract_attr(attrs, "product"),
            }
            records.append(record)

    annotations = pd.DataFrame(records)
    # Remove columns that are entirely empty
    annotations = annotations.dropna(axis=1, how="all")
    # Deduplicate by gene_id
    annotations = annotations.drop_duplicates(subset="gene_id")
    return annotations


def main():
    args = parse_args()

    print("Parsing GTF annotations...")
    annotations = parse_gtf_annotations(args.gtf)
    print(f"  Found {len(annotations)} gene annotations")

    print("Reading DESeq2 results...")
    results = pd.read_csv(args.results)
    print(f"  {len(results)} genes in results")

    print("Reading normalized counts...")
    counts = pd.read_csv(args.counts, index_col=0)
    counts.index.name = "gene_id"
    counts = counts.reset_index()
    print(f"  {len(counts)} genes in normalized counts")

    # Merge results and counts
    # The results CSV has a 'gene' column with gene IDs
    # The counts CSV has gene IDs as the row index
    gene_col = "gene" if "gene" in results.columns else results.columns[0]
    compiled = pd.merge(
        results, counts, left_on=gene_col, right_on="gene_id", how="left"
    )
    # Drop duplicate gene_id column if merge created one
    if "gene_id" in compiled.columns and gene_col != "gene_id":
        compiled = compiled.drop(columns=["gene_id"])

    # Join annotations
    compiled = pd.merge(
        annotations, compiled, left_on="gene_id", right_on=gene_col, how="right"
    )
    # Drop duplicate gene column if merge created one
    if gene_col in compiled.columns and gene_col != "gene_id":
        compiled = compiled.drop(columns=[gene_col])

    print(f"Writing compiled results ({len(compiled)} genes)...")
    compiled.to_csv(args.output, index=False)
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
