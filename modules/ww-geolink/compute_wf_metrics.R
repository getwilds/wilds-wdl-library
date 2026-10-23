#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

option_list <- list(
  make_option(dest="inputs", type="character", help="Comma-separated list of daily extract CSVs"),
  make_option(dest="pctile", type="numeric", default=0.90, help="P90 threshold"),
  make_option(dest="min_days", type="integer", default=2, help="Min days for smoke wave"),
  make_option(dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Helper: Count smoke waves (Logic from original script)
count_smoke_waves <- function(date, exposure, threshold, min_days) {
  x <- data.frame(date = as.Date(date), exposure = as.numeric(exposure)) %>%
    dplyr::arrange(date)
  if (nrow(x) == 0 || is.na(threshold)) return(NA_integer_)
  
  x$high_day <- !is.na(x$exposure) & x$exposure > threshold
  x$new_run <- c(TRUE, diff(x$date) != 1 | x$high_day[-1] != x$high_day[-nrow(x)])
  x$run_id <- cumsum(x$new_run)
  
  x %>%
    dplyr::group_by(run_id, high_day) %>%
    dplyr::summarise(run_length = dplyr::n(), .groups = "drop") %>%
    dplyr::summarise(n_smoke_waves = sum(high_day & run_length >= min_days)) %>%
    dplyr::pull(n_smoke_waves)
}

# Main function to calc annual metrics
calc_annual_metrics <- function(data, pctile, min_days) {
  data <- data %>%
    mutate(date = as.Date(date), year = year(date)) %>%
    arrange(patientID, date)

  # Weekly mean
  weekly <- data %>%
    group_by(patientID, year, isoweek(date)) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  # 1. Weeks mean > 5
  w_gt5 <- weekly %>%
    group_by(patientID, year) %>%
    summarise(n_weeks_gt5 = sum(mean_val > 5, na.rm = TRUE), .groups = "drop")

  # 2. Days > 0
  d_gt0 <- data %>%
    group_by(patientID, year) %>%
    summarise(n_days_gt0 = sum(value > 0, na.rm = TRUE), .groups = "drop")

  # 3. Peak week mean
  peak_w <- weekly %>%
    group_by(patientID, year) %>%
    summarise(peak_week_mean = max(mean_val, na.rm = TRUE), .groups = "drop")

  # 4. Smoke waves
  yearly_p90 <- data %>%
    group_by(year) %>%
    summarise(p90 = quantile(value, probs = pctile, na.rm = TRUE), .groups = "drop")

  swaves <- data %>%
    left_join(yearly_p90, by = "year") %>%
    group_by(patientID, year) %>%
    group_modify(~ {
      p90_val <- first(.x$p90)
      tibble(n_smoke_waves = count_smoke_waves(.x$date, .x$value, p90_val, min_days))
    }) %>%
    ungroup()

  # Annual lag (3 year)
  annual_mean <- data %>%
    group_by(patientID, year) %>%
    summarise(mean_daily = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(patientID, year) %>%
    group_by(patientID) %>%
    mutate(mean_3ylag = rowMeans(cbind(lag(mean_daily, 1), lag(mean_daily, 2), lag(mean_daily, 3)), na.rm = FALSE))

  # Merge all
  annual_mean %>%
    left_join(w_gt5, by = c("patientID", "year")) %>%
    left_join(d_gt0, by = c("patientID", "year")) %>%
    left_join(peak_w, by = c("patientID", "year")) %>%
    left_join(swaves, by = c("patientID", "year"))
}

# Process inputs
input_files <- unlist(strsplit(opt$inputs, ","))
daily_df <- map_df(input_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  # Assume filenames contain the date (e.g., "extract_2020-01-01.csv")
  d$date <- as.Date(gsub("extract_|.csv", "", basename(f)))
  return(d)
}) %>%
  rename(value = value)

final_metrics <- calc_annual_metrics(daily_df, opt$pctile, opt$min_days)
write_csv(final_metrics, opt$out)
