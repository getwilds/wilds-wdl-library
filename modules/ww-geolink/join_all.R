#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

option_list <- list(
  make_option("--base", dest="base", type="character", help="Base patient CSV"),
  make_option("--extracts", dest="extracts", type="character", help="Comma-separated list of extract CSVs"),
  make_option("--annuals", dest="annuals", type="character", help="Annual metrics CSV"),
  make_option("--out", dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 1. Load base
base_df <- read_csv(opt$base, show_col_types = FALSE)

# 2. Load extracts (patientID, raster, buffer_m, value) and pivot to one column per raster and buffer
extract_files <- trimws(unlist(strsplit(opt$extracts, ",")))
extracts_wide <- map_df(extract_files, function(f) read_csv(f, show_col_types = FALSE)) %>%
  mutate(col_name = paste0(raster, "_buf", buffer_m, "m")) %>%
  pivot_wider(id_cols = patientID, names_from = col_name, values_from = value)

# 3. Load annual metrics (one row per patient, year, and buffer) and pivot to one row per patient
annuals_wide <- read_csv(opt$annuals, show_col_types = FALSE) %>%
  pivot_wider(
    id_cols = patientID,
    names_from = c(year, buffer_m),
    values_from = -c(patientID, year, buffer_m),
    names_glue = "{.value}_{year}_buf{buffer_m}m"
  )

# 4. Final join
final <- base_df %>%
  left_join(extracts_wide, by = "patientID") %>%
  left_join(annuals_wide, by = "patientID")

write_csv(final, opt$out)
