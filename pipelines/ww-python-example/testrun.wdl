version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-python-example/pipelines/ww-python-example/ww-python-example.wdl" as python_example

workflow run_python_script_example {
    meta {
        author: "Taylor Firman"
        email: "tfirman@fredhutch.org"
        description: "Test workflow that generates inline test data and runs the Python script example"
        outputs: {
            results: "Output file produced by the example Python script"
        }
    }

    # Generate test data inline
    call create_test_data

    # Run the main workflow
    call python_example.run_python_script {
        input:
            python_script = create_test_data.test_script,
            sample_name = "test_sample",
            num_iterations = 3,
            input_file = create_test_data.test_input
    }

    output {
        File results = run_python_script.results
    }
}

task create_test_data {
    command <<<
        # Create a small test input file
        cat > test_input.txt <<'INPUT'
hello from the test input file
INPUT

        # Create the example Python script
        cat > test_script.py <<'SCRIPT'
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--iterations", required=True, type=int)
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.input_file) as f:
        input_data = f.read().strip()

    lines = [
        f"Sample: {args.name}",
        f"Input data: {input_data}",
        f"Iterations: {args.iterations}",
    ]
    for i in range(1, args.iterations + 1):
        lines.append(f"  iteration {i}: processed")

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Done! Wrote results to {args.output}")

if __name__ == "__main__":
    main()
SCRIPT
    >>>

    output {
        File test_script = "test_script.py"
        File test_input = "test_input.txt"
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: 1
        memory: "1 GB"
    }
}
