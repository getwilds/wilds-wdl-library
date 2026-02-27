"""Example Python script to demonstrate WDL integration.

This script accepts command-line arguments for a sample name, iteration
count, and input file, then writes a summary to an output file.
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Example script for WDL demo")
    parser.add_argument("--name", required=True, help="Sample name")
    parser.add_argument("--iterations", required=True, type=int, help="Number of iterations")
    parser.add_argument("--input-file", required=True, help="Path to input file")
    parser.add_argument("--output", required=True, help="Path to output file")
    args = parser.parse_args()

    # Read the input file
    with open(args.input_file) as f:
        input_data = f.read().strip()

    # Do some "work" with the inputs
    lines = []
    lines.append(f"Sample: {args.name}")
    lines.append(f"Input data: {input_data}")
    lines.append(f"Iterations: {args.iterations}")
    for i in range(1, args.iterations + 1):
        lines.append(f"  iteration {i}: processed")

    # Write results
    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Done! Wrote results to {args.output}")


if __name__ == "__main__":
    main()
