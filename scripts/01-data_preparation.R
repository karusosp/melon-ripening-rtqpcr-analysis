# ------------------------------- DATA PREPARATION ---------------------------
# data-preparation.R
# Combines raw RT-qPCR data (RFU and CQ) across all cultivars into two
# tidy, analysis-ready datasets saved as:
#   - cq_combined.csv
#   - rfu_combined.csv  <- (annotated with target and sample from cq_combined)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})


#-- Step 1. Loading Cq Data ---------------------------------------------------


# read_cq() is a function that takes file path of the cq dataset per cultivar
# and put the cultivar information into all the rows (well) based on the 
# name of the dataset.
# Output: cultivar | well | target | content | sample |
#         cultivar_code | stage | bio_rep | cq | set_point

read_cq <- function(path) {
  cultivar <- basename(path) |> str_remove("^cq-raw_") |> str_remove("\\.csv$")
  read_csv(path, show_col_types = FALSE,
           col_select = c(Well, Target, Content, Sample, Cq, `Set Point`)) |>
  rename(
    well      = Well,
    target    = Target,
    content   = Content,
    sample    = Sample,
    cq        = Cq,
    set_point = `Set Point`
  ) |>
    select(well, target, content, sample, cq, set_point) |>
    mutate(cultivar = cultivar)   # adding cultivar metadata 
}

# combine all the dataset for all cultivars into one table
message("reading raw cq data from 'data/cq/' directory")
cq_combined <- list.files("data/cq/", full.names = TRUE) |>
  map(read_cq) |>
  list_rbind() |>
  # Parse sample code: [CultivarPrefix].[Stage].[BiologicalRep]
  # e.g. "G.15.A" -> cultivar_code = "G", stage = "15", bio_rep = "A"
  #      "G.M.B"  -> cultivar_code = "G", stage = "M",  bio_rep = "B"
  # NTC wells have no sample code; parsed columns will be NA.
  separate(
    col    = sample,
    into   = c("cultivar_code", "stage", "bio_rep"),
    sep    = "\\.",
    remove = FALSE,
    fill   = "right"
  ) |>
  mutate(
    cultivar  = factor(cultivar),
    target    = factor(target),
    content   = factor(content),
    stage     = factor(stage, levels = c("M", "15", "30")),
    bio_rep   = factor(bio_rep,   levels = c("A", "B", "C"))
  ) |>
  select(cultivar, well, target, content, sample,
         cultivar_code, stage, bio_rep, cq, set_point)


# -- Step 2. Loading the RFU Data ---------------------------------------------
# Provide annotation to RFU dataset, primarily about sample name and gene target
# Input: Each rfu-raw_{cultivar}.csv: 40 rows (cycles) x 96 well columns (wide format)
# Output: cultivar | well | target | sample | cycle | rfu

# Well-level annotation: one row per cultivar-well (no cycle dimension)
well_annotation <- cq_combined |>
  select(cultivar, well, target, sample)

read_rfu <- function(path) {
  cultivar_name <- basename(path) |> str_remove("^rfu-raw_") |> str_remove("\\.csv$")
  read_csv(path, show_col_types = FALSE) |>
    pivot_longer(
      cols      = -Cycle,
      names_to  = "well",
      values_to = "rfu"
    ) |>
    rename(cycle = Cycle) |>
    mutate(
      well     = str_replace(well, "^([A-H])(\\d)$", "\\10\\2"), # standardize well code format so that it is aligned with cq_combined$well
      cultivar = cultivar_name #adding cultivar data based on the file name
    ) |>
    inner_join( #annotate the raw combined rfu data with more sample information for easier anlysis
      well_annotation |> filter(cultivar == cultivar_name),
      by = c("cultivar", "well")
    ) |>
    select(cultivar, well, target, sample, cycle, rfu)
}

message("reading raw rfu data from 'data/rfu/' directory")
rfu_combined <- list.files("data/rfu/", full.names = TRUE) |>
  map(read_rfu) |>
  list_rbind()

#-- Step 3. Loading Morphometric Data
morpho_df <- read_csv("data/morphometrics/fruit_morphometrics.csv", 
                      col_select = seq(1,10)) |> 
  mutate(
    samples = factor(samples)
  ) |> 
  filter(samples != is.na(samples)) |> 
  separate(col = samples,
           into = c("cultivar_code", "stage", "biorep"),
           sep = "\\.",
           remove = TRUE,
           fill = "right") |> 
  mutate(
    cultivar = case_when(
      cultivar_code == "M" ~ "Melona",
      cultivar_code == "G" ~ "GMP",
      cultivar_code == "K" ~ "Kinaya",
      cultivar_code == "T" ~ "Tacapa G"
    )
  ) |> 
  select(!cultivar_code) |> 
  relocate(cultivar)
  

#-- Step 4. Saving The Data -----------------------------------------------------------

dir.create("data/processed_data", recursive = TRUE, showWarnings = FALSE)

# Preparing the directory for analysis results
dir.create("results/")
dir.create("results/tables/") 
dir.create("results/figures/")

# write all the combined data
write_csv(cq_combined,  "data/processed_data/cq_combined.csv")
write_csv(rfu_combined,  "data/processed_data/rfu_combined.csv")
write_csv(morpho_df, "data/processed_data/morphometrics.csv")
