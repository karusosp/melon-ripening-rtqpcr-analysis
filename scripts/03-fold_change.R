# =============================================================================
# RT-qPCR Relative Gene Expression Analysis
# Method: Livak & Schmittgen (2001) — 2^(-ΔΔCq) method
# =============================================================================
#
# Experimental design:
#   - 4 melon cultivars: GMP, Kinaya, Melona, Tacapa Gold
#   - Genes of interest: CmACS, CmACO1, CmATH, CmEREBP
#   - Reference gene: CmCYP
#   - Stages: 15 DAA (calibrator/control), 30 DAA, Mature (M)
#   - 3 biological replicates per cultivar × stage (A, B, C)
#   - 2 technical replicates per well (averaged before ΔCq)
#
# Algorithm:
#   1. Average Cq of technical replicates (same gene × sample × well pairing)
#   2. ΔCq = mean_Cq(GOI) − mean_Cq(CmCYP)  [per bio-rep]
#   3. ΔΔCq = ΔCq(sample) − mean_ΔCq(calibrator, 15 DAA same cultivar)
#   4. Relative expression = 2^(−ΔΔCq)
#   5. Summary: mean ± SE across bio-reps per cultivar × stage × gene
#
# Reference:
#   Livak KJ, Schmittgen TD. Methods 2001;25:402–8.
#   doi:10.1006/meth.2001.1262

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

# -- Parameters ---------------------------------------------------------------

REFERENCE_GENE <- "CmCYP"
CALIBRATOR     <- "15"   # stage used as reference (15 DAA)
GOI            <- c("CmACS", "CmACO", "CmATH", "CmEREBP")

# Cultivar prefix → full name mapping
CULTIVAR_MAP <- c(
  G = "GMP",
  K = "Kinaya",
  M = "Melona",
  T = "Tacapa"
)

# -- Step 1. Load and combine raw Cq files ------------------------------------

read_cq_raw <- function(path) {
  cultivar_code <- basename(path) |>
    str_remove("cq-raw_") |>
    str_remove("\\.csv$")

  read_csv(path, show_col_types = FALSE) |>
    rename_with(str_trim) |>                    # strip any BOM-induced spaces
    rename(
      well    = Well,
      target  = Target,
      content = Content,
      sample  = Sample,
      cq      = Cq
    ) |>
    select(well, target, content, sample, cq) |>
    filter(content != "NTC") |>  # drop NTC and failed wells
    separate(
      col    = sample,
      into   = c("cultivar_prefix", "stage", "bio_rep"),
      sep    = "\\.",
      remove = FALSE
    ) |>
    mutate(
      cultivar = CULTIVAR_MAP[cultivar_prefix],
      stage    = factor(stage, levels = c("15", "30", "M")),
      bio_rep  = factor(bio_rep, levels = c("A", "B", "C"))
    )
}

message("Loading raw Cq files ...")
cq_files <- list.files("data/cq/",
                        full.names = TRUE)

cq_all <- map(cq_files, read_cq_raw) |> 
  list_rbind()

message(sprintf("  Loaded %d Cq measurements across %d cultivars",
                nrow(cq_all), n_distinct(cq_all$cultivar)))


# -- Step 2. Average technical replicates -------------------------------------
# Each bio-rep × gene combination may appear in two wells (technical replicates).
# We average them here before any ΔCq calculation.

cq_mean <- cq_all |>
  group_by(cultivar, stage, bio_rep, target) |>
  summarise(
    mean_cq     = mean(cq, na.rm = TRUE),
    n_tech_reps = n(),
    .groups     = "drop"
  )


# -- Step 3. Pivot to wide and compute ΔCq ------------------------------------
# Wide format: one row per cultivar × stage × bio_rep,
# columns = one mean_Cq per gene

cq_wide <- cq_mean |>
  pivot_wider(
    names_from  = target,
    values_from = mean_cq
  )

# ΔCq = Cq(GOI) − Cq(reference)
# Computed per bio-rep; 

delta_cq <- cq_wide |>
  mutate(across(
    all_of(GOI),
    ~ . - .data[[REFERENCE_GENE]],
    .names = "dCq_{.col}"
  )) |>
  select(cultivar, stage, bio_rep, starts_with("dCq_"))


# -- Step 4. Extract calibrator ΔCq per cultivar × bio_rep × gene -------------
 
calibrator_dCq <- delta_cq |>
  filter(stage == CALIBRATOR) |>
  select(cultivar, bio_rep, starts_with("dCq_")) |>
  rename_with(~ sub("^dCq_", "calib_dCq_", .x), starts_with("dCq_"))
 
# Join on cultivar + bio_rep so each sample row gets its paired 15 DAA ΔCq
delta_cq <- delta_cq |>
  left_join(calibrator_dCq, by = c("cultivar", "bio_rep"))
  

# -- Step 5. ΔΔCq and fold-change (2^−ΔΔCq) per bio-rep -----------------------

# Build ΔΔCq columns dynamically for each GOI
for (gene in GOI) {
  dcq_col   <- paste0("dCq_",    gene)
  calib_col <- paste0("calib_dCq_", gene)
  ddcq_col  <- paste0("ddCq_",   gene)
  fc_col    <- paste0("fc_",     gene)

  delta_cq[[ddcq_col]] <- delta_cq[[dcq_col]] - delta_cq[[calib_col]]
  delta_cq[[fc_col]]   <- 2^(-delta_cq[[ddcq_col]])
}

# Clean up helper columns
results_per_biorep <- delta_cq |>
  select(cultivar, stage, bio_rep,
         starts_with("dCq_"),
         starts_with("ddCq_"),
         starts_with("fc_"))


# -- Step 6. Summary table: mean ± SE per cultivar × stage × gene -------------

fc_cols <- paste0("fc_", GOI)

results_summary <- results_per_biorep |>
  pivot_longer(
    cols      = all_of(fc_cols),
    names_to  = "gene",
    values_to = "fold_change"
  ) |>
  mutate(gene = str_remove(gene, "^fc_")) |>
  group_by(cultivar, stage, gene) |>
  summarise(
    n_bio_reps  = sum(!is.na(fold_change)),
    mean_fc     = mean(fold_change, na.rm = TRUE),
    se_fc       = sd(fold_change, na.rm = TRUE) / sqrt(n_bio_reps),
    sd_fc       = sd(fold_change, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(
    gene  = factor(gene,  levels = GOI),
    stage = factor(stage, levels = c("15", "30", "M")),
    cultivar = factor(cultivar, levels = c("Melona", "GMP", "Kinaya", "Tacapa"))
  ) |>
  arrange(cultivar, gene, stage)

# Also keep a long-form per-biorep table (useful for stats / plotting)
results_long <- results_per_biorep |>
  pivot_longer(
    cols      = all_of(fc_cols),
    names_to  = "gene",
    values_to = "fold_change"
  ) |>
  mutate(gene = str_remove(gene, "^fc_")) |>
  select(cultivar, stage, bio_rep, gene, fold_change) |>
  mutate(
    gene     = factor(gene,     levels = GOI),
    stage    = factor(stage,    levels = c("15", "30", "M")),
    cultivar = factor(cultivar, levels = c("Melona", "GMP", "Kinaya", "Tacapa"))
  ) |>
  arrange(cultivar, gene, stage, bio_rep)

# -- Step 7. Saving the Fold Change Results -----------------------------------
write_csv(results_long,    "results/tables/ddcq_per_biorep.csv")
write_csv(results_summary, "results/tables/ddcq_summary.csv")

