# =============================================================================
# RT-qPCR Efficiency Estimation
# Method: Linear regression on log10(RFU) per cycle — Ramakers et al. (2003)
# =============================================================================
#
# Algorithm (per well):
#   1. Log10-transform RFU (Relative Fluoresence Unit)
#   2. Slide a window of `window_size` cycles across all available cycles
#   3. For each window, fit lm(log10_RFU ~ Cycle) and record slope and R²
#   4. Discard windows with slope < min_slope to remove noise region
#   5. Among eligible windows, pick the one with the highest R²
#   6. Efficiency per well: E = 10^best_slope
#   7. Summarize mean efficiency per target gene across all wells
#
# The calculation for one well at a time is conducted by
# `estimate_efficiency_one_well` function.
# And that function run reiteratively for the whole 96-well dataset 
# by `estimate_efficiency_full` function.

# Reference:
#   Ramakers et al. (2003) Neuroscience Letters 339:62-66
#   doi:10.1016/S0304-3940(02)01423-4


suppressWarnings(
  library(tidyverse) 
)

# ── Parameters ────────────────────────────────────────────────────────────────

window_size <- 5        # cycles per window (Ramakers recommends 4–6)
minimal_slope   <- 0.15 # minimum slope to exclude plateau/noise regions
                        # log10 scale: 0.15 corresponds to E >= ~1.41
cultivars   <- c("gmp", "kinaya", "melona", "tacapa")


# ── Core Function: estimate_efficiency_one_well() ───────────────────────────────────
# The function find the best linear window and calculate its efficiency per well.
# It takes a data frame and for a single well (columns: Cycle, log10_RFU),
# returns the best slope, R², and window boundaries for that well.

estimate_efficiency_one_well <- function(df_well, window = window_size,
                             min_slope = minimal_slope) {

  cycles <- df_well$cycle
  vals   <- df_well$log10_rfu
  n      <- length(cycles)

  # --- Step 1: fit lm for every possible window position ---
  all_windows <- map(seq_len(n - window + 1), function(i) {
    idx <- i:(i + window - 1)
    fit <- lm(vals[idx] ~ cycles[idx])
    tibble(
      cycle_start = cycles[idx[1]],
      cycle_end   = cycles[idx[window]],
      slope       = coef(fit)[[2]],
      r2          = summary(fit)$r.squared
    )
  }) |>
    bind_rows()

  # --- Step 2: exclude plateau and noise regions ---
  # A minimum slope of 0.15 (log10 scale) corresponds to E >= ~1.41
  # This removes flat plateau windows 
  eligible <- filter(all_windows, slope >= min_slope)

  # Guard: if no windows pass the filter, return NA
  if (nrow(eligible) == 0) {
    return(tibble(slope = NA_real_, efficiency = NA_real_,
                  r2 = NA_real_, cycle_start = NA_real_, cycle_end = NA_real_))
  }

  # --- Step 3: pick the window with the best R² ---
  # Best R² = most linear region (The best candidate for exponential phase)
  best <- slice_max(eligible, order_by = r2, n = 1)

  # --- Step 4: calculate efficiency from slope ---
  # E = 10^slope because we used log10 (per Ramakers et al.)
  mutate(best, efficiency = 10^slope)
}


# ── Wrapper: Estimate Efficiency for a Full RFU dataset ──────────────────

# Takes the joined/annotated RFU + Cq data frame 
# Returns one row per well with efficiency metrics

estimate_efficiency_full <- function(df) {
  df |>
    filter(rfu > 0,
           sample != "NA") |> # filter out the NTCs
    mutate(log10_rfu = log10(rfu)) |>
    group_by(cultivar, well, target, sample) |>
    group_modify(~ estimate_efficiency_one_well(.x)) |>
    ungroup()
}


# ── Summary: Mean Efficiency per Target Gene ──────────────────────────────────
#
# Collapses per-well results to one row per gene
# Flags primers deviating from the ideal 90–110% efficiency range (E = 1.9–2.1)

summarize_efficiency <- function(df) {
  df |>
    group_by(target) |>
    summarise(
      n_wells         = n(),
      mean_efficiency = mean(efficiency, na.rm = TRUE),
      sd_efficiency   = sd(efficiency,   na.rm = TRUE),
      mean_r2         = mean(r2,         na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      efficiency_pct  = round((mean_efficiency - 1) * 100, 1),
      efficiency_flag = case_when(
        mean_efficiency >= 1.90 & mean_efficiency <= 2.10 ~ "OK",
        mean_efficiency >= 1.80 & mean_efficiency <  1.90 ~ "Acceptable (low)",
        mean_efficiency >  2.10 & mean_efficiency <= 2.20 ~ "Acceptable (high)",
        TRUE                                              ~ "Check primer"
      )
    )
}


# ── Load data ─────────────────────────────────────────────────────────────────

RFU_PATH <- "data/processed_data/rfu_combined.csv"
CQ_PATH  <- "data/processed_data/cq_combined.csv"

message("Reading combined data files ...")
rfu_combined <- read_csv(RFU_PATH) 
cq_combined  <- read_csv(CQ_PATH) 

# ── Run ───────────────────────────────────────────────────────────────────────

# One row per well per cultivar with slope, R2, cycle window, and efficiency
message("Estimating efficiency per well ...")
efficiency_per_well <- estimate_efficiency_full(rfu_combined)

# Per-gene summary collapsed across all cultivars and wells
message("Summarizing efficiency per gene targets")
efficiency_per_gene <- summarize_efficiency(efficiency_per_well)

# Per-gene per-cultivar: collapse per-well results within each cultivar
message("Summarizing efficiency per gene per cultivar ...")
efficiency_per_gene_per_cultivar <- efficiency_per_well |>
  group_by(cultivar, target) |>
  summarise(
    n_wells         = n(),
    mean_efficiency = mean(efficiency, na.rm = TRUE),
    sd_efficiency   = sd(efficiency,   na.rm = TRUE),
    mean_r2         = mean(r2,         na.rm = TRUE),
    .groups         = "drop"
  ) |>
  mutate(
    efficiency_pct  = round((mean_efficiency - 1) * 100, 1),
    efficiency_flag = case_when(
      mean_efficiency >= 1.90 & mean_efficiency <= 2.10 ~ "OK",
      mean_efficiency >= 1.80 & mean_efficiency <  1.90 ~ "Acceptable (low)",
      mean_efficiency >  2.10 & mean_efficiency <= 2.20 ~ "Acceptable (high)",
      TRUE                                              ~ "Check primer"
    )
  ) |>
  arrange(cultivar, target)


# ── Save ──────────────────────────────────────────────────────────────────────

dir.create("results", showWarnings = FALSE)

write_csv(efficiency_per_well,              "results/tables/efficiency-per-well.csv")
write_csv(efficiency_per_gene,              "results/tables/efficiency-per-gene.csv")
write_csv(efficiency_per_gene_per_cultivar, "results/tables/efficiency-per-gene-per-cultivar.csv")

message("Done.")
message(sprintf("  efficiency-per-gene.csv              : %d genes",    nrow(efficiency_per_gene)))
message(sprintf("  efficiency-per-well.csv              : %d wells",    nrow(efficiency_per_well)))
message(sprintf("  efficiency-per-gene-per-cultivar.csv : %d rows",     nrow(efficiency_per_gene_per_cultivar)))

