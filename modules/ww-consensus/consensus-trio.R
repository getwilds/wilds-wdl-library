#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("4 arguments must be provided: GATK file, SAMtools file, Mutect file, and base file name.", call.=FALSE)
}

GATKfile <- args[1]
SAMfile <- args[2]
MutectFile <- args[3]
baseName <- args[4]

library(tidyverse)

# Simple function to read Annovar-annotated table and extract key variant info
read_variants <- function(filepath, caller_name) {
  df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE)

  # Create variant ID from position columns (these should be consistent across Annovar output)
  df$VariantID <- paste(df$Chr, df$Start, df$End, df$Ref, df$Alt, sep = "-")
  df$Caller <- caller_name

  return(df)
}

print("Reading GATK variants...")
gatk <- read_variants(GATKfile, "GATK")
print(paste("  Found", nrow(gatk), "variants"))

print("Reading SAMtools variants...")
sam <- read_variants(SAMfile, "SAMtools")
print(paste("  Found", nrow(sam), "variants"))

print("Reading Mutect variants...")
mutect <- read_variants(MutectFile, "Mutect")
print(paste("  Found", nrow(mutect), "variants"))

# Get unique variant IDs from each caller
gatk_ids <- unique(gatk$VariantID)
sam_ids <- unique(sam$VariantID)
mutect_ids <- unique(mutect$VariantID)

# Find all unique variants across all callers
all_ids <- unique(c(gatk_ids, sam_ids, mutect_ids))
print(paste("Total unique variants across all callers:", length(all_ids)))

# Build consensus table
consensus <- data.frame(VariantID = all_ids, stringsAsFactors = FALSE)

# Add caller detection flags
consensus$InGATK <- consensus$VariantID %in% gatk_ids
consensus$InSAMtools <- consensus$VariantID %in% sam_ids
consensus$InMutect <- consensus$VariantID %in% mutect_ids

# Count how many callers detected each variant
consensus$CallerCount <- rowSums(consensus[, c("InGATK", "InSAMtools", "InMutect")])

# Assign confidence tiers based on caller agreement
consensus$Confidence <- case_when(
  consensus$CallerCount == 3 ~ "HighConfidence",
  consensus$CallerCount == 2 ~ "MediumConfidence",
  consensus$CallerCount == 1 ~ "LowConfidence"
)

# Split VariantID back into components
consensus <- consensus %>%
  separate(VariantID, into = c("Chr", "Start", "End", "Ref", "Alt"),
           sep = "-", remove = FALSE)

# Determine variant type
consensus$Type <- ifelse(nchar(consensus$Ref) == nchar(consensus$Alt), "SNV", "INDEL")

# Add annotation data from first available caller (prefer GATK > SAMtools > Mutect)
# Only include columns that exist in the source data
get_annotations <- function(variant_id, gatk_df, sam_df, mutect_df) {
  # Try GATK first
  gatk_row <- gatk_df[gatk_df$VariantID == variant_id, ]
  if (nrow(gatk_row) > 0) return(gatk_row[1, ])

  # Then SAMtools
  sam_row <- sam_df[sam_df$VariantID == variant_id, ]
  if (nrow(sam_row) > 0) return(sam_row[1, ])

  # Finally Mutect
  mutect_row <- mutect_df[mutect_df$VariantID == variant_id, ]
  if (nrow(mutect_row) > 0) return(mutect_row[1, ])

  return(NULL)
}

# Get common annotation columns (Annovar standard columns)
annovar_cols <- c("Func.refGene", "Gene.refGene", "GeneDetail.refGene",
                  "ExonicFunc.refGene", "AAChange.refGene")

# Find which annotation columns actually exist
available_cols <- annovar_cols[annovar_cols %in% colnames(gatk)]

if (length(available_cols) > 0) {
  # Create annotation lookup from GATK (or fall back to others)
  all_variants <- bind_rows(
    gatk %>% select(VariantID, all_of(available_cols)),
    sam %>% select(VariantID, all_of(available_cols)),
    mutect %>% select(VariantID, all_of(available_cols))
  ) %>%
    distinct(VariantID, .keep_all = TRUE)

  consensus <- left_join(consensus, all_variants, by = "VariantID")
}

# Reorder columns for cleaner output
output_cols <- c("VariantID", "Chr", "Start", "End", "Ref", "Alt", "Type",
                 "InGATK", "InSAMtools", "InMutect", "CallerCount", "Confidence")
if (length(available_cols) > 0) {
  output_cols <- c(output_cols, available_cols)
}
consensus <- consensus %>% select(all_of(output_cols))

# Summary stats
print("")
print("=== Consensus Summary ===")
print(paste("High confidence (3 callers):", sum(consensus$Confidence == "HighConfidence")))
print(paste("Medium confidence (2 callers):", sum(consensus$Confidence == "MediumConfidence")))
print(paste("Low confidence (1 caller):", sum(consensus$Confidence == "LowConfidence")))
print(paste("Total variants:", nrow(consensus)))

# Write output
output_file <- paste0(baseName, ".consensus.tsv")
write.table(consensus, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Output written to:", output_file))
