#=======================================================================================
# Gene Expression Data Statistical Test
# ======================================================================================
# Determine whether or not there is significance different in gene expression across developmental stage 
# in each cultivar and each of the target gene. 
# 


#-- Step 1. ddcq statistical test
library(broom) # a package for make tidy dataframe from stats object 

tukey_results <- read_csv("results/tables/ddcq_per_biorep.csv") |>
  group_by(cultivar, gene) |>
  group_modify(function(df, keys) {
    # Drop unused factor levels so checks and contrasts reflect observed data only
    df <- df |>
      filter(!is.na(fold_change))
    # log2 transform fold-change for ANOVA (stabilizes variance)
    df <- df |> 
      mutate(log2_fc = log2(fold_change + 1e-6))
    # anova and tukey's post-hoc
    aov_fit <- aov(log2_fc ~ stage, data = df)
    tukey   <- TukeyHSD(aov_fit)

    as_tibble(tukey$stage, rownames = "comparison") |>
      rename(estimate = diff, conf.low = lwr, conf.high = upr,
             adj.p.value = `p adj`)
  }) |>
  ungroup() |>
  mutate(
    sig = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01  ~ "**",
      adj.p.value < 0.05  ~ "*",
      TRUE                ~ "ns"
    )
  )   



# -- Step 2. Saving Results
write_csv(tukey_results,   "results/tables/tukey_anova.csv")
message("  results/tables/tukey_anova.csv       — Tukey HSD pairwise stage comparisons")