version 1.0

workflow run_python_script {
    meta {
        author: "Taylor Firman"
        email: "tfirman@fredhutch.org"
        description: "Example WDL workflow that wraps a Python script, demonstrating how to pass string, integer, and file inputs as command-line arguments"
    }

    parameter_meta {
        python_script: "The Python script to execute"
        sample_name: "A string argument passed to the script via --name"
        num_iterations: "An integer argument passed to the script via --iterations"
        input_file: "A file argument passed to the script via --input-file"
    }

    input {
        File python_script
        String sample_name
        Int num_iterations
        File input_file
    }

    call run_script {
        input:
            python_script = python_script,
            sample_name = sample_name,
            num_iterations = num_iterations,
            input_file = input_file
    }

    output {
        File results = run_script.results
    }
}

task run_script {
    input {
        File python_script
        String sample_name
        Int num_iterations
        File input_file
    }

    command <<<
        python3 ~{python_script} \
            --name ~{sample_name} \
            --iterations ~{num_iterations} \
            --input-file ~{input_file} \
            --output results.txt
    >>>

    output {
        File results = "results.txt"
    }

    runtime {
        docker: "getwilds/flax:0.1.0"
        cpu: 1
        memory: "2 GB"
    }
}
