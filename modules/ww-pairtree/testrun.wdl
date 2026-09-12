version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "./ww-pairtree.wdl" as ww_pairtree
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow pairtree_example {
  # Auto-generate a small synthetic SSM + params.json for testing purposes
  call ww_testdata.create_pairtree_data { }

  # Cluster mutations into subclones
  call ww_pairtree.cluster_variants { input:
      ssm_file = create_pairtree_data.ssm_file,
      params_file = create_pairtree_data.params_file
  }

  # Sample clone trees consistent with the clustered data
  call ww_pairtree.run_pairtree { input:
      ssm_file = create_pairtree_data.ssm_file,
      params_file = cluster_variants.clustered_params_file,
      trees_per_chain = 100,
      parallel_chains = 1
  }

  # Generate an interactive visualization of the results
  call ww_pairtree.plot_tree { input:
      ssm_file = create_pairtree_data.ssm_file,
      params_file = cluster_variants.clustered_params_file,
      results_file = run_pairtree.results_file
  }

  output {
    File clustered_params_file = cluster_variants.clustered_params_file
    File results_file = run_pairtree.results_file
    File tree_html = plot_tree.tree_html
    File tree_json = plot_tree.tree_json
  }
}
