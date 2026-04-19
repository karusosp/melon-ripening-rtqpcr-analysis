# -------------------- Gene Expression Data Visualization ---------------------
# Nested Bar Chart showing the relative expression of each gene in th Y-axis
# between different developmental stage for each cultivar (shown in X-axis). 
# The Bar Chart is then annotated with CLD (compact letter display) to show the
# significance between groups (developmental stage).  

library(tidyverse)
library(multcompView)
library(patchwork)

# -- Step 1. Data Loading -----------------------------------------------------

# Fold Change Data
tidy_fc <- read_csv("results/tables/ddcq_per_biorep.csv") |> 
  mutate(
    cultivar = str_to_lower(cultivar),
    log2_fc = log2(fold_change),
    stage    = factor(stage, levels = c("15", "30", "M"))
  )

# -- Step 2: Derive CLD from Tukey pairwise results ---------------------------
library(multcompView)

cld_letter <- tidy_fc |>
  group_by(cultivar, gene) |>
  group_modify(function(df, keys) {
    aov_fit     <- aov(log2_fc ~ stage, data = df)         
    tukey       <- TukeyHSD(aov_fit)
    p_vec       <- tukey$stage[, "p adj"]
    cld_letters <- multcompLetters(p_vec)$Letters
    tibble(stage = names(cld_letters), cld_letters = cld_letters) 
  }) |>
  ungroup()

# -- Step 3: Build the summary df and join CLD labels -------------------------

plot_df <- tidy_fc |>
  group_by(cultivar, stage, gene) |>
  summarize(
    mean_fc = mean(log2_fc),
    sem     = sd(log2_fc)/sqrt(3),
    .groups = "drop"
  ) |>
  left_join(cld_letter, by = c("cultivar", "stage", "gene")) |> 
  mutate(cultivar = case_when(
    cultivar == "gmp" ~ "GMP",
    cultivar == "tacapa" ~ "Tacapa",
    cultivar == "kinaya" ~ "Kinaya",
    cultivar == "melona" ~ "Melona"
  ))

# -- Step 4: Plot with geom_text annotation -----------------------------------
log2_fc_plot <- plot_df |>
  mutate(
    cultivar = factor(cultivar, levels = c("Melona", "GMP", "Tacapa", "Kinaya")),
    gene  = factor(gene, levels = c("CmACS", "CmACO1", "CmATH", "CmEREBP")
    )
  ) |> 
  ggplot(aes(
    x     = cultivar,
    y     = mean_fc,
    group = stage,
    fill  = cultivar
  )) +
  facet_wrap(~ gene, scales = "free") +
  geom_col(position = "dodge", color = "black", width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_fc, ymax = mean_fc + sem),
    position = position_dodge(0.8),
    width    = 0.2,
    color    = "darkred"
  ) +
  # CLD letters above error bars
  geom_text(
    aes(
      y     = mean_fc + sem + 0.30,
      label = cld_letters,
      group = stage
    ),
    position = position_dodge(0.8),
    size     = 3.5,
    vjust    = 0
  ) +
  # Stage labels below the x-axis
  geom_text(
    aes(
      y     = -Inf,
      label = stage,
      group = stage
    ),
    position = position_dodge(0.8),
    size     = 3,
    vjust    = 1,          # pushes text below the axis line
    color    = "black",
    fontface = "bold"
  ) +
  scale_fill_viridis_d()+
  # scale_fill_brewer(palette = "Set1") +
  coord_cartesian(clip = "off") +   # allows text to render outside plot area
  theme_minimal(base_family = "serif") +
  theme(
    axis.title.y = element_text(colour = "black", face = "bold"),
    legend.title = element_blank(),
    text         = element_text(colour = "black"),
    strip.text   = element_text(colour = "black", face = "bold"), # bold the text per plot within the facet_wrap()
    plot.margin  = margin(t = 5, r = 5, b = 30, l = 5),   # extra bottom margin for stage labels
    axis.text.x = element_blank(),
    panel.spacing.y = unit(0.5, "cm")     # add spacing per plot panel
  ) +
  ylab("Relative Expression Log2 (Fold Change)") +
  xlab("") -> log2_fc_plot

# -- Saving The Plot ----------------------------------------------------------
ggsave(
  plot = log2_fc_plot,
  filename = "results/figures/log2_fold-change.png",
  create.dir = TRUE,
  dpi = 300,
  units = "mm",
  width = 174,
  height = 140)

