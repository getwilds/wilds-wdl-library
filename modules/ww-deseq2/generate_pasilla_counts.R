#!/usr/bin/env Rscript

# Generate individual STAR-format count files using the pasilla Bioconductor dataset
# Author: WILDS Team
# Description: Creates individual ReadsPerGene.out.tab files for each sample to mimic STAR output

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(pasilla)
    library(DESeq2)
})

# Define command line options
option_list <- list(
    make_option(c("-n", "--nsamples"), type = "integer", default = 7,
                help = "Number of samples to include (default: 7, max: 7 for pasilla dataset)"),
    make_option(c("-g", "--ngenes"), type = "integer", default = 10000,
                help = "Approximate number of genes to include (will select top expressed genes)"),
    make_option(c("-c", "--condition"), type = "character", default = "condition",
                help = "Name for the condition column in metadata"),
    make_option(c("-p", "--prefix"), type = "character", default = "pasilla",
                help = "Prefix for output files")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                          description = "Generate individual STAR-format count files from pasilla dataset")
opt <- parse_args(opt_parser)

# Main execution
tryCatch({
    cat("=== Generating individual STAR-format count files from pasilla data ===\n")
    cat("Parameters:\n")
    cat("  Samples:", opt$nsamples, "\n")
    cat("  Genes:", opt$ngenes, "\n")
    cat("  Condition column:", opt$condition, "\n")
    cat("  Output prefix:", opt$prefix, "\n\n")
    
    # Get the actual pasilla data files
    count_file <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla")
    sample_file <- system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla")
    
    cat("Loading pasilla data from:\n")
    cat("Count file:", count_file, "\n")
    cat("Sample file:", sample_file, "\n\n")
    
    # Read the count matrix and sample data
    count_data <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
    sample_data <- read.csv(sample_file, header = TRUE, row.names = 1)
    
    cat("Original data dimensions:\n")
    cat("Genes:", nrow(count_data), "\n")
    cat("Samples:", ncol(count_data), "\n\n")
    
    # Remove 'fb' suffix from sample metadata row names to match count data
    rownames(sample_data) <- gsub("fb$", "", rownames(sample_data))

    # Make sure sample names match between count data and metadata
    common_samples <- intersect(colnames(count_data), rownames(sample_data))
    count_data <- count_data[, common_samples]
    sample_data <- sample_data[common_samples, ]
    
    # Limit to requested number of samples
    n_samples_requested <- opt$nsamples
    if (n_samples_requested > ncol(count_data)) {
        n_samples_requested <- ncol(count_data)
        cat("Warning: Requested", opt$nsamples, "samples, but only", ncol(count_data), "available\n")
    }
    
    # Select samples (ensure we have both conditions represented)
    untreated_samples <- rownames(sample_data)[sample_data$condition == "untreated"]
    treated_samples <- rownames(sample_data)[sample_data$condition == "treated"]
    
    # Take roughly equal numbers from each condition
    n_untreated <- ceiling(n_samples_requested / 2)
    n_treated <- n_samples_requested - n_untreated
    
    # Adjust if we don't have enough samples of one type
    if (n_untreated > length(untreated_samples)) {
        n_untreated <- length(untreated_samples)
        n_treated <- n_samples_requested - n_untreated
    }
    if (n_treated > length(treated_samples)) {
        n_treated <- length(treated_samples)
        n_untreated <- n_samples_requested - n_treated
    }
    
    selected_samples <- c(
        untreated_samples[1:n_untreated],
        treated_samples[1:n_treated]
    )
    
    # Subset data
    count_data <- count_data[, selected_samples]
    sample_data <- sample_data[selected_samples, ]
    
    # Select top expressed genes
    n_genes_requested <- opt$ngenes
    if (n_genes_requested > nrow(count_data)) {
        n_genes_requested <- nrow(count_data)
        cat("Warning: Requested", opt$ngenes, "genes, but only", nrow(count_data), "available\n")
    }
    
    # Calculate mean counts per gene and select top expressed
    gene_means <- rowMeans(count_data)
    top_genes <- order(gene_means, decreasing = TRUE)[1:n_genes_requested]
    
    count_data <- count_data[top_genes, ]
    
    # Create individual STAR-format count files for each sample
    individual_count_files <- character(0)
    sample_names <- character(0)
    sample_conditions <- character(0)
    
    for (i in 1:ncol(count_data)) {
        sample_name <- colnames(count_data)[i]
        sample_condition <- as.character(sample_data[sample_name, "condition"])
        
        # Generate filename in STAR format: samplename.ReadsPerGene.out.tab
        count_filename <- paste0(sample_name, ".ReadsPerGene.out.tab")
        
        # Get counts for this sample
        sample_counts <- count_data[, i]
        
        # Calculate some realistic summary statistics for STAR header
        # These are fake but representative values
        total_reads <- sum(sample_counts) * 2  # Approximate total reads
        uniquely_mapped <- round(sum(sample_counts) * 0.85)  # ~85% uniquely mapped
        multimapping <- round(sum(sample_counts) * 0.10)     # ~10% multimapping
        unmapped <- total_reads - uniquely_mapped - multimapping
        
        # Create STAR-format output with header statistics
        # STAR ReadsPerGene.out.tab format:
        # First 4 lines are summary statistics
        # Then: gene_id, unstranded_count, stranded_forward, stranded_reverse
        
        # Open file for writing
        file_conn <- file(count_filename, "w")
        
        # Write STAR header lines (summary statistics)
        writeLines(paste("N_unmapped", unmapped, sep = "\t"), file_conn)
        writeLines(paste("N_multimapping", multimapping, sep = "\t"), file_conn)
        writeLines(paste("N_noFeature", round(total_reads * 0.02), sep = "\t"), file_conn)  # ~2% no feature
        writeLines(paste("N_ambiguous", round(total_reads * 0.03), sep = "\t"), file_conn)  # ~3% ambiguous
        
        # Write gene counts
        # Format: gene_id, unstranded, forward_strand, reverse_strand
        # We'll use the actual counts as "unstranded" and generate reasonable strand-specific counts
        for (j in 1:length(sample_counts)) {
            gene_id <- rownames(count_data)[j]
            unstranded_count <- sample_counts[j]
            
            # Generate strand-specific counts (simulate roughly 60/40 split for strand specificity)
            forward_count <- round(unstranded_count * runif(1, 0.3, 0.7))
            reverse_count <- unstranded_count - forward_count
            
            writeLines(paste(gene_id, unstranded_count, forward_count, reverse_count, sep = "\t"), file_conn)
        }
        
        close(file_conn)
        
        # Store information for outputs
        individual_count_files <- c(individual_count_files, count_filename)
        sample_names <- c(sample_names, sample_name)
        sample_conditions <- c(sample_conditions, sample_condition)
        
        cat("Created:", count_filename, "for sample", sample_name, "with condition", sample_condition, "\n")
    }
    
    # Create sample names file
    sample_names_file <- paste0(opt$prefix, "_sample_names.txt")
    writeLines(sample_names, sample_names_file)
    
    # Create sample conditions file
    sample_conditions_file <- paste0(opt$prefix, "_sample_conditions.txt")
    writeLines(sample_conditions, sample_conditions_file)
    
    # Create count files list
    count_files_list <- paste0(opt$prefix, "_count_files.txt")
    writeLines(individual_count_files, count_files_list)
    
    # Write gene information
    gene_info_file <- paste0(opt$prefix, "_gene_info.txt")
    gene_info <- data.frame(
        gene_id = rownames(count_data),
        gene_name = rownames(count_data),  # Use gene ID as name since we don't have symbols
        stringsAsFactors = FALSE
    )
    write.table(gene_info, 
                file = gene_info_file, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE,
                quote = FALSE)
    
    # Print summary
    cat("\n=== Generated individual STAR-format count files ===\n")
    cat("Files created:\n")
    for (i in 1:length(individual_count_files)) {
        cat("  Count file:", individual_count_files[i], "(", sample_names[i], "-", sample_conditions[i], ")\n")
    }
    cat("  Sample names list:", sample_names_file, "\n")
    cat("  Sample conditions list:", sample_conditions_file, "\n")
    cat("  Count files list:", count_files_list, "\n")
    cat("  Gene info:", gene_info_file, "\n\n")
    cat("Data summary:\n")
    cat("  Total samples:", length(sample_names), "\n")
    cat("  Genes per file:", nrow(count_data), "\n")
    cat("  Conditions:", paste(unique(sample_conditions), collapse = ", "), "\n")
    cat("  Sample names:", paste(sample_names, collapse = ", "), "\n")
    cat("\nScript completed successfully!\n")

}, error = function(e) {
    cat("Error occurred during execution:\n")
    cat(conditionMessage(e), "\n")
    quit(status = 1)
})
