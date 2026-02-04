#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("4 arguments must be provided: HaplotypeCaller file, mpileup file, Mutect file, and base file name.", call.=FALSE)
}

HaploFile <- args[1]
MpileupFile <- args[2]
MutectFile <- args[3]
baseName <- args[4]

library(tidyverse)

# Simple function to read Annovar-annotated table and extract key variant info
# Note: Annovar's Otherinfo column contains extra tab-separated VCF fields
# beyond the defined header columns, so we use read_tsv which handles this gracefully
read_variants <- function(filepath, caller_name) {
  df <- read_tsv(filepath, show_col_types = FALSE, col_types = cols(.default = "c"))

  # Create variant ID from position columns (these should be consistent across Annovar output)
  df$VariantID <- paste(df$Chr, df$Start, df$End, df$Ref, df$Alt, sep = "-")
  df$Caller <- caller_name

  return(df)
}

print("Reading HaplotypeCaller variants...")
haplo <- read_variants(HaploFile, "HaplotypeCaller")
print(paste("  Found", nrow(haplo), "variants"))

print("Reading mpileup variants...")
mpileup <- read_variants(MpileupFile, "mpileup")
print(paste("  Found", nrow(mpileup), "variants"))

print("Reading Mutect variants...")
mutect <- read_variants(MutectFile, "Mutect")
print(paste("  Found", nrow(mutect), "variants"))

# Get unique variant IDs from each caller
haplo_ids <- unique(haplo$VariantID)
mpileup_ids <- unique(mpileup$VariantID)
mutect_ids <- unique(mutect$VariantID)

# Find all unique variants across all callers
all_ids <- unique(c(haplo_ids, mpileup_ids, mutect_ids))
print(paste("Total unique variants across all callers:", length(all_ids)))

# Build consensus table
consensus <- data.frame(VariantID = all_ids, stringsAsFactors = FALSE)

# Add caller detection flags
consensus$InHaplotypeCaller <- consensus$VariantID %in% haplo_ids
consensus$InMpileup <- consensus$VariantID %in% mpileup_ids
consensus$InMutect <- consensus$VariantID %in% mutect_ids

# Count how many callers detected each variant
consensus$CallerCount <- rowSums(consensus[, c("InHaplotypeCaller", "InMpileup", "InMutect")])

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

# Add annotation data from first available caller (prefer HaplotypeCaller > mpileup > Mutect)
# Only include columns that exist in the source data
get_annotations <- function(variant_id, haplo_df, mpileup_df, mutect_df) {
  # Try HaplotypeCaller first
  haplo_row <- haplo_df[haplo_df$VariantID == variant_id, ]
  if (nrow(haplo_row) > 0) return(haplo_row[1, ])

  # Then mpileup
  mpileup_row <- mpileup_df[mpileup_df$VariantID == variant_id, ]
  if (nrow(mpileup_row) > 0) return(mpileup_row[1, ])

  # Finally Mutect
  mutect_row <- mutect_df[mutect_df$VariantID == variant_id, ]
  if (nrow(mutect_row) > 0) return(mutect_row[1, ])

  return(NULL)
}

# Get common annotation columns (Annovar standard columns)
annovar_cols <- c("Func.refGene", "Gene.refGene", "GeneDetail.refGene",
                  "ExonicFunc.refGene", "AAChange.refGene")

# Find which annotation columns actually exist
available_cols <- annovar_cols[annovar_cols %in% colnames(haplo)]

if (length(available_cols) > 0) {
  # Create annotation lookup from HaplotypeCaller (or fall back to others)
  all_variants <- bind_rows(
    haplo %>% select(VariantID, all_of(available_cols)),
    mpileup %>% select(VariantID, all_of(available_cols)),
    mutect %>% select(VariantID, all_of(available_cols))
  ) %>%
    distinct(VariantID, .keep_all = TRUE)

  consensus <- left_join(consensus, all_variants, by = "VariantID")
}

# Reorder columns for cleaner output
output_cols <- c("VariantID", "Chr", "Start", "End", "Ref", "Alt", "Type",
                 "InHaplotypeCaller", "InMpileup", "InMutect", "CallerCount", "Confidence")
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
