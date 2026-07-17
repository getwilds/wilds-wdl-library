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


ATTR_KEYS = ["gene_name", "gene", "gene_biotype", "product", "locus_tag", "db_xref"]


def parse_gtf_annotations(gtf_path):
    """Parse gene annotations from a GTF file.

    Walks every feature row (gene/transcript/exon/CDS/etc.), keys by gene_id,
    and takes the first non-null value seen per attribute. Handles GTFs that
    lack `gene` rows (e.g. UCSC ncbiRefSeq) and NCBI's `gene "name"` attribute
    (used in place of Ensembl's `gene_name`).
    """
    by_gene = {}
    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            gene_id = extract_attr(attrs, "gene_id")
            if not gene_id:
                continue
            if gene_id not in by_gene:
                by_gene[gene_id] = {"gene_id": gene_id}
            for key in ATTR_KEYS:
                if by_gene[gene_id].get(key):
                    continue
                value = extract_attr(attrs, key)
                if value:
                    by_gene[gene_id][key] = value

    if not by_gene:
        return pd.DataFrame(columns=["gene_id"])

    # NCBI GTFs use `gene "name"` instead of `gene_name "name"`; promote it.
    for gene_id in by_gene:
        if not by_gene[gene_id].get("gene_name") and by_gene[gene_id].get("gene"):
            by_gene[gene_id]["gene_name"] = by_gene[gene_id]["gene"]

    annotations = pd.DataFrame(list(by_gene.values()))
    annotations = annotations.drop(columns=["gene"], errors="ignore")
    # Remove columns that are entirely empty, but always keep gene_id
    keep_cols = [c for c in annotations.columns
                 if c == "gene_id" or annotations[c].notna().any()]
    annotations = annotations[keep_cols]
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
    # Standardize the gene identifier column to 'gene_id'
    if gene_col != "gene_id":
        compiled = compiled.rename(columns={gene_col: "gene_id"})
    # Drop any duplicate gene_id columns from the counts side
    compiled = compiled.loc[:, ~compiled.columns.duplicated()]

    # Join annotations
    compiled = pd.merge(
        annotations, compiled, on="gene_id", how="right"
    )

    print(f"Writing compiled results ({len(compiled)} genes)...")
    compiled.to_csv(args.output, index=False)
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
