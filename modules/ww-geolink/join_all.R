#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(readr)

option_list <- list(
  make_option(dest="base", type="character", help="Base patient CSV"),
  make_option(dest="extracts", type="character", help="Comma-separated list of extract CSVs"),
  make_option(dest="annuals", type="character", help="Annual metrics CSV"),
  make_option(dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load base
base_df <- read_csv(opt$base, show_col_types = FALSE)

# 2. Load and pivot extracts
# We expect extract CSVs to be: patientID, value
extract_files <- unlist(strsplit(opt$extracts, ","))
all_ext_df <- map_df(extract_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  col_name <- gsub("extract_|.csv", "", basename(f))
  d %>% rename(!!col_name := value)
}) %>%
  # Since extract_val.R outputs patientID and value, we have multiple patientID cols
  # We need to join them by patientID
  reduce(left_join, by = "patientID")

# 3. Load annuals
annuals_df <- read_csv(opt$annuals, show_col_types = FALSE)

# 4. Final Join
final <- base_df %>%
  left_join(all_ext_df, by = "patientID") %>%
  left_join(annuals_df, by = "patientID")

write_csv(final, opt$out)
