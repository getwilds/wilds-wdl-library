#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(lubridate)
library(readr)
library(purrr)
library(tidyr)

option_list <- list(
  make_option("--inputs", dest="inputs", type="character", help="Comma-separated list of daily extract CSVs"),
  make_option("--pctile", dest="pctile", type="numeric", default=0.90, help="Percentile for smoke-wave threshold"),
  make_option("--min_days", dest="min_days", type="integer", default=2, help="Min days for smoke wave"),
  make_option("--out", dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Helper: count smoke waves (runs of at least min_days consecutive days above threshold)
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

# Calculate annual metrics for a single buffer distance
calc_annual_metrics <- function(data, pctile, min_days) {
  data <- data %>%
    mutate(year = year(date), week = isoweek(date)) %>%
    arrange(patientID, date)

  # Weekly mean
  weekly <- data %>%
    group_by(patientID, year, week) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  # 1. Weeks with mean > 5
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

  # 5. Annual mean with 3-year lag
  annual_mean <- data %>%
    group_by(patientID, year) %>%
    summarise(mean_daily = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(patientID, year) %>%
    group_by(patientID) %>%
    mutate(mean_3ylag = rowMeans(cbind(lag(mean_daily, 1), lag(mean_daily, 2), lag(mean_daily, 3)), na.rm = FALSE)) %>%
    ungroup()

  # Merge all
  annual_mean %>%
    left_join(w_gt5, by = c("patientID", "year")) %>%
    left_join(d_gt0, by = c("patientID", "year")) %>%
    left_join(peak_w, by = c("patientID", "year")) %>%
    left_join(swaves, by = c("patientID", "year"))
}

# Load all extracts; the date is parsed from the raster name (YYYY-MM-DD or YYYYMMDD)
input_files <- trimws(unlist(strsplit(opt$inputs, ",")))
daily_df <- map_df(input_files, function(f) read_csv(f, show_col_types = FALSE))

date_str <- regmatches(daily_df$raster, regexpr("[0-9]{4}-?[0-9]{2}-?[0-9]{2}", daily_df$raster))
if (length(date_str) != nrow(daily_df)) {
  stop("Could not parse a date (YYYY-MM-DD or YYYYMMDD) from every raster name")
}
daily_df$date <- as.Date(gsub("-", "", date_str), format = "%Y%m%d")

# Metrics are computed separately for each buffer distance
final_metrics <- daily_df %>%
  group_split(buffer_m) %>%
  map_df(function(d) {
    res <- calc_annual_metrics(d, opt$pctile, opt$min_days)
    res$buffer_m <- d$buffer_m[1]
    res
  })

write_csv(final_metrics, opt$out)
