## STIC immune surveillance: 2x2 factorial mixed-effects model (IFN tertile x Metabolic tertile)
##


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
})

base_dir <- "D:/CyCIF_Geomx_analysis_LSP_HMS/CYCIF_ORION_atlas_2.0/Xenium_analysis"
in_csv   <- file.path(base_dir, "immune_survellance_score", "STIC_immune_surveillance_score_joint_IFN_Metab.csv")
out_dir  <- file.path(base_dir, "reproduce_mixed_effects_model")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metrics      <- c("Density_z", "Function_z")
joint_levels <- c("LowIFN_LowMetab", "LowIFN_HighMetab", "HighIFN_LowMetab", "HighIFN_HighMetab")

## 1. Load and derive the two 2-level factors ---------------------------------
pm <- read.csv(in_csv, stringsAsFactors = FALSE)
stopifnot(all(pm$joint_group %in% joint_levels), all(metrics %in% names(pm)))

pm <- pm %>%
  mutate(
    joint_group   = factor(joint_group, levels = joint_levels),
    IFN_tertile   = factor(ifelse(grepl("^HighIFN", joint_group), "High", "Low"), levels = c("Low", "High")),
    Metab_tertile = factor(ifelse(grepl("HighMetab$", joint_group), "High", "Low"), levels = c("Low", "High"))
  )

cat(sprintf("Input: %d rows, %d patients\n", nrow(pm), n_distinct(pm$sample_id)))
print(table(pm$IFN_tertile, pm$Metab_tertile, dnn = c("IFN", "Metab")))
if (any(pm$n_cells < 5)) {
  cat("\nLow-confidence patient x group cells (n_cells < 5):\n")
  print(pm %>% filter(n_cells < 5) %>% select(sample_id, joint_group, n_cells))
}

## 2. Mixed models: treatment coding (simple effects) and sum coding (main effects)
fit_2x2 <- function(metric, coding = c("treatment", "sum")) {
  coding <- match.arg(coding)
  d <- pm
  if (coding == "sum") {
    contrasts(d$IFN_tertile)   <- contr.sum(2)
    contrasts(d$Metab_tertile) <- contr.sum(2)
  } else {
    contrasts(d$IFN_tertile)   <- contr.treatment(2)
    contrasts(d$Metab_tertile) <- contr.treatment(2)
  }
  m  <- lmer(as.formula(paste0(metric, " ~ IFN_tertile * Metab_tertile + (1 | sample_id)")), data = d)
  cf <- coef(summary(m))
  term_names <- c("IFN", "Metab", "IFN x Metab")
  effect_type <- if (coding == "treatment")
    c("simple effect: IFN High vs Low within Metab-Low",
      "simple effect: Metab High vs Low within IFN-Low",
      "interaction") else
    c("main effect (averaged over Metab); estimate = half the High-Low difference, sign: Low minus High",
      "main effect (averaged over IFN); estimate = half the High-Low difference, sign: Low minus High",
      "interaction")
  data.frame(
    metric      = metric,
    coding      = coding,
    term        = term_names,
    effect_type = effect_type,
    estimate    = unname(cf[2:4, "Estimate"]),
    std_error   = unname(cf[2:4, "Std. Error"]),
    df          = unname(cf[2:4, "df"]),
    t_value     = unname(cf[2:4, "t value"]),
    p_value     = unname(cf[2:4, "Pr(>|t|)"]),
    singular_fit = isSingular(m),
    stringsAsFactors = FALSE
  )
}

model_res <- bind_rows(lapply(metrics, function(mt) bind_rows(fit_2x2(mt, "treatment"), fit_2x2(mt, "sum"))))

cat("\n================ 2x2 mixed model: metric ~ IFN * Metab + (1|sample_id) ================\n")
print(model_res %>% select(metric, coding, term, estimate, std_error, df, p_value, singular_fit) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 4))), row.names = FALSE)
write.csv(model_res, file.path(out_dir, "2x2_mixed_model_results.csv"), row.names = FALSE)

## 3. Group means +/- SEM ------------------------------------------------------
group_stats <- pm %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  group_by(metric, joint_group, IFN_tertile, Metab_tertile) %>%
  summarise(n_patients = sum(!is.na(value)),
            mean = mean(value, na.rm = TRUE),
            sem  = sd(value, na.rm = TRUE) / sqrt(n_patients),
            .groups = "drop") %>%
  arrange(metric, joint_group)
cat("\nGroup means +/- SEM:\n")
print(group_stats %>% mutate(across(c(mean, sem), ~ round(.x, 3))), n = Inf)
write.csv(group_stats, file.path(out_dir, "2x2_group_means_SEM.csv"), row.names = FALSE)

## 4. Per-patient paired Wilcoxon contrasts ------------------------------------
wide <- pm %>%
  select(sample_id, joint_group, all_of(metrics)) %>%
  pivot_wider(names_from = joint_group, values_from = all_of(metrics), names_sep = "__")

contrast_defs <- list(
  "IFN effect @ Metab Low"               = c("HighIFN_LowMetab",  "LowIFN_LowMetab"),
  "IFN effect @ Metab High"              = c("HighIFN_HighMetab", "LowIFN_HighMetab"),
  "Metab effect @ IFN Low"               = c("LowIFN_HighMetab",  "LowIFN_LowMetab"),
  "Metab effect @ IFN High"              = c("HighIFN_HighMetab", "HighIFN_LowMetab"),
  "Combined (double-High vs double-Low)" = c("HighIFN_HighMetab", "LowIFN_LowMetab")
)

paired_res <- bind_rows(lapply(names(contrast_defs), function(cn) {
  bind_rows(lapply(metrics, function(mt) {
    hi <- wide[[paste0(mt, "__", contrast_defs[[cn]][1])]]
    lo <- wide[[paste0(mt, "__", contrast_defs[[cn]][2])]]
    ok <- !is.na(hi) & !is.na(lo)
    dif <- hi[ok] - lo[ok]
    data.frame(contrast = cn, metric = mt, n_patients = sum(ok),
               mean_diff = mean(dif), n_increase = sum(dif > 0), n_decrease = sum(dif < 0), n_tie = sum(dif == 0),
               p_value = wilcox.test(hi[ok], lo[ok], paired = TRUE)$p.value,
               stringsAsFactors = FALSE)
  }))
}))
cat("\nPer-patient paired Wilcoxon contrasts:\n")
print(paired_res %>% mutate(across(c(mean_diff, p_value), ~ signif(.x, 4))), row.names = FALSE)
write.csv(paired_res, file.path(out_dir, "2x2_paired_wilcoxon_contrasts.csv"), row.names = FALSE)


