library(tidyverse)
library(readxl)
library(seriation)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(hrbrthemes)
library(viridis)
library(powerjoin)
library(lme4)
library(emmeans)

# ── Output directory ──────────────────────────────────────────────────────────
out_dir <- "CD8_subset_quantification"
dir.create(file.path(out_dir, "panckpos"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "panckneg"), recursive = TRUE, showWarnings = FALSE)

# ── Populations ───────────────────────────────────────────────────────────────
# Note: leading space trimmed from "CD8apCD103n"
numerator_markers <- c(
  "CD8apPD1pKi67pLAG3nTIM3n",
  "CD8apLAG3pTIM3n",
  "CD8apLAG3pTIM3p",
  "CD8apCD103n",
  "CD8apTbetpNKG2ApCD103p",
  "CD8apTbetpNKG2ApCD103n",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p"
)

pop_labels <- c(
  "CD8apPD1pKi67pLAG3nTIM3n"                   = "CD8a+PD1+Ki67+LAG3-TIM3- (prolif, non-exhausted) / total CD8a+ (%)",
  "CD8apLAG3pTIM3n"                             = "CD8a+LAG3+TIM3- / total CD8a+ (%)",
  "CD8apLAG3pTIM3p"                             = "CD8a+LAG3+TIM3+ (exhausted) / total CD8a+ (%)",
  "CD8apCD103n"                                 = "CD8a+CD103- (non-resident) / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103p"                      = "CD8a+Tbet+NKG2A+CD103+ / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103n"                      = "CD8a+Tbet+NKG2A+CD103- / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n"  = "CD8a+Tbet+NKG2A+CD103+PD1+Ki67+LAG3-TIM3- / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n"  = "CD8a+Tbet+NKG2A+CD103-PD1+Ki67+LAG3-TIM3- / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103pLAG3pTIM3p"           = "CD8a+Tbet+NKG2A+CD103+LAG3+TIM3+ / total CD8a+ (%)",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p"           = "CD8a+Tbet+NKG2A+CD103-LAG3+TIM3+ / total CD8a+ (%)"
)

denominator_marker <- "CD8ap"

# ── GLMM with observational-level random effect (OLRE) ───────────────────────
# OLRE (1|obs_id) accounts for overdispersion in binomial models by giving
# each observation its own random intercept, absorbing extra-binomial variation.
# (1|slideName) accounts for multiple ROIs from the same patient/slide.
run_glmm <- function(data, pop_name) {
  df <- data %>%
    mutate(
      obs_id = as.character(row_number()),
      n_neg  = pmax(n_denominator - n_numerator, 0)
    ) %>%
    filter(!is.na(n_numerator), !is.na(n_denominator), n_denominator > 0)

  tryCatch({
    model <- glmer(
      cbind(n_numerator, n_neg) ~ ROIcat + (1 | slideName) + (1 | obs_id),
      family  = binomial,
      data    = df,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    model_null <- glmer(
      cbind(n_numerator, n_neg) ~ 1 + (1 | slideName) + (1 | obs_id),
      family  = binomial,
      data    = df,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    lrt   <- anova(model_null, model)
    lrt_p <- lrt$`Pr(>Chisq)`[2]

    emm <- emmeans(model, ~ ROIcat)
    pwc <- as.data.frame(pairs(emm, adjust = "BH")) %>%
      mutate(
        population = pop_name,
        sig = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE             ~ "ns"
        )
      )

    list(model = model, lrt_p = lrt_p, pairs = pwc)

  }, error = function(e) {
    warning("GLMM failed for ", pop_name, ": ", e$message)
    NULL
  })
}

# ── Summary statistics helper ─────────────────────────────────────────────────
compute_summary_stats <- function(data) {
  data %>%
    group_by(population_numerator, ROIcat) %>%
    summarise(
      n       = n(),
      min     = min(proportion,    na.rm = TRUE),
      max     = max(proportion,    na.rm = TRUE),
      median  = median(proportion, na.rm = TRUE),
      mean    = mean(proportion,   na.rm = TRUE),
      sd      = sd(proportion,     na.rm = TRUE),
      .groups = "drop"
    )
}

# ── Strip plot helper ─────────────────────────────────────────────────────────
make_strip_plot <- function(data_pop, pop, label, plot_title = "epithelial-immune subset") {
  glmm <- run_glmm(data_pop, pop)

  subtitle_text <- if (!is.null(glmm)) {
    sprintf("GLMM (OLRE+slide RE): ROIcat p = %.3g", glmm$lrt_p)
  } else {
    "GLMM: failed to converge"
  }

  y_max <- max(data_pop$proportion, na.rm = TRUE) * 1.15

  p <- ggplot(data_pop,
              aes(x = ROIcat, y = proportion, color = Incidental_numerator)) +
    geom_jitter(width = 0.2, height = 0, size = 4, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 color = "black", linewidth = 0.3) +
    scale_color_manual(
      values = c("Y" = "black", "N" = "red"),
      name   = "Incidental"
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, y_max),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    labs(
      y        = label,
      x        = "morphology",
      title    = plot_title,
      subtitle = subtitle_text
    ) +
    theme_classic() +
    theme(
      axis.text     = element_text(color = "black"),
      axis.title    = element_text(color = "black"),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  list(plot = p, glmm = glmm)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: panckpos (PanCK+, epithelial compartment)
# ══════════════════════════════════════════════════════════════════════════════
marker_data <- read_csv("panckposmean.csv")

marker_data_long <- marker_data %>%
  select(
    slideName, ROIname_n, ROIcat, ROIname, cell_count,
    BRCA, Incidental,
    matches("^mean_[A-Za-z0-9_-]+[pn]$")
  ) %>%
  pivot_longer(
    cols      = matches("^mean_[A-Za-z0-9_-]+[pn]$", ignore.case = FALSE),
    names_to  = "population",
    values_to = "proportion"
  ) %>%
  filter(!population %in% c("mean_ORF1p", "mean_Orientation")) %>%
  mutate(
    population = str_replace(population, fixed("mean_"), ""),
    n          = round(cell_count * proportion)
  )

marker_data_long_subset <- marker_data_long %>%
  filter(!ROIcat %in% c("", "NA"))

marker_data_long_subset$ROIcat <- factor(
  marker_data_long_subset$ROIcat,
  levels = c("normalFT", "normalFimbria", "p53signature", "STIL", "STIC", "cancer")
)

marker_data_numerators  <- marker_data_long_subset %>% filter(population %in% numerator_markers)
marker_data_denominator <- marker_data_long_subset %>% filter(population == denominator_marker)

marker_data_proportion_pos <- power_inner_join(
  marker_data_numerators,
  marker_data_denominator,
  by     = c("slideName", "ROIname", "ROIname_n", "ROIcat"),
  check  = check_specs(
    unmatched_keys_left  = "warn",
    unmatched_keys_right = "warn",
    duplicate_keys_right = "warn"
  ),
  suffix = c("_numerator", "_denominator")
) %>%
  mutate(proportion = n_numerator / n_denominator)

# Summary statistics — panckpos
summary_stats_pos <- compute_summary_stats(marker_data_proportion_pos)
write_csv(summary_stats_pos, file.path(out_dir, "panckpos", "summary_stats_panckpos.csv"))

# Loop over all populations — panckpos
glmm_results_pos <- list()

for (pop in numerator_markers) {
  data_pop <- marker_data_proportion_pos %>% filter(population_numerator == pop)
  result   <- make_strip_plot(data_pop, pop, pop_labels[[pop]],
                              plot_title = "panckpos: epithelial-immune subset")

  ggsave(
    filename = file.path(out_dir, "panckpos", paste0(pop, "_panckpos.pdf")),
    plot     = result$plot,
    width    = 8, height = 5
  )
  ggsave(
    filename = file.path(out_dir, "panckpos", paste0(pop, "_panckpos.png")),
    plot     = result$plot,
    width    = 8, height = 5, dpi = 300
  )

  if (!is.null(result$glmm)) glmm_results_pos[[pop]] <- result$glmm$pairs
  message("panckpos done: ", pop)
}

# Save pairwise GLMM results — panckpos
if (length(glmm_results_pos) > 0) {
  bind_rows(glmm_results_pos) %>%
    write_csv(file.path(out_dir, "panckpos", "GLMM_pairwise_panckpos.csv"))
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: panckneg (PanCK-, stromal/immune compartment)
# ══════════════════════════════════════════════════════════════════════════════
marker_data <- read_csv("pancknegmean.csv")

marker_data_long <- marker_data %>%
  select(
    slideName, ROIname_n, ROIcat, ROIname, cell_count,
    BRCA, Incidental,
    matches("^mean_[A-Za-z0-9_-]+[pn]$")
  ) %>%
  pivot_longer(
    cols      = matches("^mean_[A-Za-z0-9_-]+[pn]$", ignore.case = FALSE),
    names_to  = "population",
    values_to = "proportion"
  ) %>%
  filter(!population %in% c("mean_ORF1p", "mean_Orientation")) %>%
  mutate(
    population = str_replace(population, fixed("mean_"), ""),
    n          = round(cell_count * proportion)
  )

marker_data_long_subset <- marker_data_long %>%
  filter(!ROIcat %in% c("", "NA"))

marker_data_long_subset$ROIcat <- factor(
  marker_data_long_subset$ROIcat,
  levels = c("normalFT", "normalFimbria", "p53signature", "STIL", "STIC", "cancer")
)

marker_data_numerators  <- marker_data_long_subset %>% filter(population %in% numerator_markers)
marker_data_denominator <- marker_data_long_subset %>% filter(population == denominator_marker)

marker_data_proportion_neg <- power_inner_join(
  marker_data_numerators,
  marker_data_denominator,
  by     = c("slideName", "ROIname", "ROIname_n", "ROIcat"),
  check  = check_specs(
    unmatched_keys_left  = "warn",
    unmatched_keys_right = "warn",
    duplicate_keys_right = "warn"
  ),
  suffix = c("_numerator", "_denominator")
) %>%
  mutate(proportion = n_numerator / n_denominator)

# Summary statistics — panckneg
summary_stats_neg <- compute_summary_stats(marker_data_proportion_neg)
write_csv(summary_stats_neg, file.path(out_dir, "panckneg", "summary_stats_panckneg.csv"))

# Loop over all populations — panckneg
glmm_results_neg <- list()

for (pop in numerator_markers) {
  data_pop <- marker_data_proportion_neg %>% filter(population_numerator == pop)
  result   <- make_strip_plot(data_pop, pop, pop_labels[[pop]],
                              plot_title = "panckneg: stromal-immune subset")

  ggsave(
    filename = file.path(out_dir, "panckneg", paste0(pop, "_panckneg.pdf")),
    plot     = result$plot,
    width    = 8, height = 5
  )
  ggsave(
    filename = file.path(out_dir, "panckneg", paste0(pop, "_panckneg.png")),
    plot     = result$plot,
    width    = 8, height = 5, dpi = 300
  )

  if (!is.null(result$glmm)) glmm_results_neg[[pop]] <- result$glmm$pairs
  message("panckneg done: ", pop)
}

# Save pairwise GLMM results — panckneg
if (length(glmm_results_neg) > 0) {
  bind_rows(glmm_results_neg) %>%
    write_csv(file.path(out_dir, "panckneg", "GLMM_pairwise_panckneg.csv"))
}

message("All done. Plots and GLMM results saved to: ", out_dir)
