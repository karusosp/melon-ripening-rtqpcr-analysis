# =============================================================================
#                         MORPHOMETRICS DATA ANALYSIS
# Statistical approach:
#   - One-way ANOVA (stage) per trait per cultivar
#   - Post-hoc Tukey HSD within each cultivar across stages
#   - Normality is assumed based on biological plausibility 
#
# Outputs:
#   results/morpho_summary.csv           — mean ± SD per trait per cultivar
#   results/morpho_anova.csv             — one-way ANOVA per trait
#   results/morpho_tukey_ssc.csv         — Tukey HSD for SSC
#   results/morpho_tukey_firmness.csv    — Tukey HSD for firmness
#   figures/morpho_ssc_firmness.pdf      — publication figure
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(multcompView)
  library(ggpubr)#for CLD annotation 
})

#-- Step 1. Data Loading ------------------------------------------------------
morpho_df <- read_csv("data/morphometrics/fruit_morphometrics.csv", 
               show_col_types = FALSE, 
               col_select = seq(1,10), 
               n_max = 36 #select only non-NA cell
               ) 

#-- Step 2. Decoding Sample Code and Tidying The Data -------------------------
morpho_df <- morpho_df |> 
  separate( # Decoding sample code inside `samples` column
    col = samples,
    into = c("cultivar", "stage", "biorep"), 
    sep = "\\.",
    remove = TRUE,
    fill = "right"
  ) |>  
  mutate(
    cultivar = case_when( # Adding proper label
      cultivar == "M" ~ "Melona", 
      cultivar == "K" ~ "Kinaya",
      cultivar == "T" ~ "Tacapa", 
      cultivar == "G" ~ "GMP"
    ),
    stage    = case_when(
      stage == "15" ~ "15 DAA", 
      stage == "30" ~ "30 DAA", 
      stage == "M"  ~ "Mature"
    )
  ) |> 
  mutate(
    cultivar = factor(cultivar, levels = c("Melona", "GMP", "Tacapa", "Kinaya")),
    stage    = factor(stage, levels = c("15 DAA", "30 DAA", "Mature"))
  )

#-- Step 2. Summary Statistics Table — mean ± SD per cultivar × stage ---------
summary_tbl <- morpho_df %>%
  group_by(cultivar, stage) %>%
  summarise(
    across(
      where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

summary_tbl <- summary_tbl %>%
  mutate(across(
    ends_with("_mean"),
    ~ paste0(
      round(.x, 2),
      " (± ",
      round(get(sub("_mean", "_sd", cur_column())), 2),
      ")"
    ),
    .names = "{.col}_fmt"
  ))

summary_tbl <- summary_tbl |> 
  select(cultivar, stage, ends_with("_mean_fmt")) |> 
  rename_with(~ gsub("_mean_fmt", "", .x)) |> 
  rename("fw" = fruit_weight, 
         "vc" = v_perimeter,
         "hc" = h_perimeter,
         "vd" = v_diameter,
         "hd" = h_diameter,
         "rd" = rind_thickness,
         "ft" = flesh_thickness,
         "ssc"= ssc,
         "firm" = firmness)

#-- Step 3. One-Way ANOVA and Tukey Test-----------------------------------------------------
morpho_anova <- morpho_df |> 
  pivot_longer(
    cols = c(fruit_weight, v_perimeter, 
             h_perimeter, v_diameter, 
             h_diameter, 
             flesh_thickness, ssc, firmness),
    names_to = "traits",
    values_to = "values"
  ) |> 
  mutate(traits = factor(traits)) |> 
  group_by(traits, cultivar) |> 
  group_modify(function(df, keys){
    aov_fit <- aov(values ~ stage, data = df)
    tukey <- TukeyHSD(aov_fit)
    as_tibble(tukey$stage, rownames = "comparison")
  }) |> 
  ungroup()

#--Step 4. Compact Letter Display for Tukey Test ------------------------------
morpho_cld <- morpho_df |> 
  pivot_longer(
    cols = c(fruit_weight, v_perimeter, 
             h_perimeter, v_diameter, 
             h_diameter, 
             flesh_thickness, ssc, firmness),
    names_to = "traits",
    values_to = "values"
  ) |> 
  mutate(traits = factor(traits)) |> 
  group_by(traits, cultivar) |> 
  group_modify(function(df, keys){
    aov_fit <- aov(values ~ stage, data = df)
    tukey <- TukeyHSD(aov_fit)
    p_vec <- tukey$stage[, "p adj"]
    cld_letters <- multcompLetters(p_vec)$Letters
    tibble(
           stage = names(cld_letters), 
           cld_letters = cld_letters)
  })
 
#-- Step 5. Firmness and SSC Plot ---------------------------------------------
traits <- c("fruit_weight", "v_perimeter", "h_perimeter",
            "v_diameter", "h_diameter",
            "flesh_thickness", "ssc", "firmness")

morpho_plot <- morpho_df |> 
  pivot_longer(
    cols = c(fruit_weight, v_perimeter, 
             h_perimeter, v_diameter, 
             h_diameter, 
             flesh_thickness, ssc, firmness),
    names_to = "traits",
    values_to = "values"
  ) |> 
  filter(
    traits %in% c("firmness", "ssc")
  ) |> 
  group_by(traits, cultivar, stage) |> 
  summarize(
    mean_value = mean(values), 
    sd_value = sd(values)
  ) |> 
  ggplot(
    aes(
      x = stage,
      y = mean_value, 
      color = cultivar,
    )
  )+
  facet_wrap(~traits, 
             scale = "free_y",
             labeller = as_labeller(c(firmness = "Firmness (N)",
                                      ssc = "Soluble Solid Content (\u00B0Brix)")))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean_value - sd_value,
        ymax = mean_value + sd_value),
    width = 0.1)+
  geom_line(aes(group = cultivar))+
  scale_color_viridis_d()+
  theme_pubclean()+
  labs(
    x = "",
    y = ""
  )+
  theme(
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 10),
  )

#-- Step 5. Saving The Data and Plot

# Saving Tables
write_csv(morpho_anova, "results/tables/morpho_anova.csv")
write_csv(morpho_cld, "results/tables/morpho_cld.csv")
write_csv(summary_tbl, "results/tables/morpho_summary.csv")

# Saving Plot
ggsave(
  plot = morpho_plot,
  filename = "results/figures/morpho_plot.png",
  create.dir = TRUE,
  dpi = 300,
  units = "mm",
  width = 174,
  height = 140
  )

