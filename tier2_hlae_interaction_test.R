# tier2_hlae_interaction_test.R
# Formal test of whether the HLA-E+ vs HLA-E- immune-neighborhood effect
# (within LDH+FASN+ index cells) DIFFERS between Normal/p53sig and STIC.
#
# Follow-up to tier2_hlae_fimbria_p53_vs_stic.R, which fit the two histology
# groups as SEPARATE models -- a descriptive comparison only (non-overlapping
# CIs, no formal test of the group difference itself). This script instead
# fits ONE joint model per immune population with an explicit
# HLA_Ep_idx x roi_group interaction, and tests that interaction directly
# via an emmeans "ratio of odds ratios" (ROR) contrast:
#   ROR = OR(HLA-E+ vs HLA-E-, in STIC) / OR(HLA-E+ vs HLA-E-, in Normal/p53sig)
# ROR > 1 (CI excluding 1) means the HLA-E effect is significantly STRONGER
# in STIC; ROR < 1 means significantly weaker (i.e. stronger in Normal/p53sig).
#
# Same data restriction and same per-slide balance filter as the parent
# script (>=5 cells in both HLA-E- and HLA-E+, applied within each
# roi_group independently) -- uses the identical 24 (Normal/p53sig) + 16
# (STIC) slide set.
#
# Model, per immune population:
#   cbind(cnt, fails) ~ HLA_Ep_idx * roi_group + (1 | slideName)
# BH-FDR applied across the 14 immune populations on the interaction p-value.

library(tidyverse)
library(arrow)
library(lme4)
library(emmeans)
library(scales)

# ── Settings ──────────────────────────────────────────────────────────────────
nbhd_parquet        <- "STIC_analysis/data_with_neighborhoods_50um_filt.parquet"
out_dir             <- "STIC_analysis/results/neighborhood_metabolic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_slide <- 5   # same threshold as tier2_hlae_fimbria_p53_vs_stic.R

immune_pops <- c(
  "CD8ap", "CD4p", "NKG2ApCD8anCD4n", "CD11cp", "CD68p",
  "CD8apLAG3pTIM3p", "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p", "CD4pTIM3pLAG3p",
  "CD11cpCD68pCD103nTIM3p", "CD11cpCD68nCD103nTIM3p",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n", "CD68pCD11cp"
)
immune_cnt_cols <- paste0(immune_pops, "_count")

roi_levels <- c("Normal / p53sig", "STIC")

# ── Load & filter (identical restriction to tier2_hlae_fimbria_p53_vs_stic.R) ─
cat("Reading neighborhood parquet ...\n")
df_raw <- read_parquet(nbhd_parquet) %>%
  rename(HLA_Ep_idx = `HLA-Ep`)
cat(sprintf("  %d cells, %d slides\n", nrow(df_raw), n_distinct(df_raw$slideName)))

df <- df_raw %>%
  filter(
    PanCKp     == TRUE,
    Incidental == "Y",
    ROIcat %in% c("normalFimbria", "p53signature", "STIC"),
    LDHp       == TRUE,
    FASNp      == TRUE,
    !is.na(HLA_Ep_idx)
  ) %>%
  mutate(
    roi_group = case_when(
      ROIcat %in% c("normalFimbria", "p53signature") ~ "Normal / p53sig",
      ROIcat == "STIC"                                ~ "STIC"
    ),
    roi_group  = factor(roi_group, levels = roi_levels),
    HLA_Ep_idx = factor(HLA_Ep_idx, levels = c(FALSE, TRUE), labels = c("HLA-E-", "HLA-E+"))
  )

slide_counts <- df %>%
  count(slideName, roi_group, HLA_Ep_idx) %>%
  pivot_wider(names_from = HLA_Ep_idx, values_from = n, values_fill = 0L)

slides_ok <- slide_counts %>%
  filter(`HLA-E-` >= min_cells_per_slide, `HLA-E+` >= min_cells_per_slide) %>%
  distinct(slideName, roi_group)

df <- df %>% semi_join(slides_ok, by = c("slideName", "roi_group"))

cat("\nCell counts after slide filter (used in interaction GLMM):\n")
print(table(df$roi_group, df$HLA_Ep_idx))
cat("Slides used per roi_group:\n")
print(df %>% group_by(roi_group) %>% summarise(n_slides = n_distinct(slideName)))

# ── Interaction GLMM helper: one immune population ──────────────────────────────
run_interaction_glmm <- function(cnt_col, data) {
  pop <- sub("_count$", "", cnt_col)
  tryCatch({
    grp <- data %>%
      select(slideName, HLA_Ep_idx, roi_group, total_neighbors, cnt = all_of(cnt_col)) %>%
      filter(!is.na(cnt)) %>%
      mutate(fails = total_neighbors - cnt)

    m <- glmer(
      cbind(cnt, fails) ~ HLA_Ep_idx * roi_group + (1 | slideName),
      family  = binomial,
      data    = grp,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )

    em <- emmeans(m, ~ HLA_Ep_idx * roi_group)

    # per-group simple HLA-E+ vs HLA-E- odds ratios, from the JOINT model
    # (for reference / sanity check against the separately-fit models)
    simple <- contrast(em, "trt.vs.ctrl", ref = 1, by = "roi_group",
                       type = "response", infer = c(TRUE, TRUE)) %>%
      as_tibble() %>%
      select(roi_group, odds.ratio, asymp.LCL, asymp.UCL) %>%
      mutate(immune_pop = pop)

    # interaction contrast: ratio of odds ratios (STIC vs Normal/p53sig)
    ror <- contrast(em, interaction = c("trt.vs.ctrl1", "trt.vs.ctrl1"),
                    type = "response", infer = c(TRUE, TRUE)) %>%
      as_tibble() %>%
      rename(ROR = odds.ratio, ROR.LCL = asymp.LCL, ROR.UCL = asymp.UCL) %>%
      mutate(immune_pop = pop)

    list(simple = simple, ror = ror)
  }, error = function(e) {
    warning("Interaction GLMM failed for ", pop, ": ", conditionMessage(e))
    NULL
  })
}

cat("\nFitting interaction GLMMs (14 pops) ...\n")
res_list <- map(immune_cnt_cols, run_interaction_glmm, data = df) %>% compact()

simple_res <- map_dfr(res_list, "simple") %>%
  mutate(immune_pop = factor(immune_pop, levels = immune_pops)) %>%
  arrange(immune_pop, roi_group)

ror_res <- map_dfr(res_list, "ror") %>%
  mutate(
    immune_pop = factor(immune_pop, levels = immune_pops),
    p_adj      = p.adjust(p.value, method = "BH"),
    sig_label  = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    significant = p_adj < 0.05,
    direction   = case_when(
      !significant       ~ "no significant difference",
      ROR > 1             ~ "significantly stronger in STIC",
      ROR < 1             ~ "significantly stronger in Normal/p53sig"
    )
  ) %>%
  arrange(p_adj)

cat("\n=== Interaction test: ratio of odds ratios (STIC vs Normal/p53sig) ===\n")
ror_res %>%
  select(immune_pop, ROR, ROR.LCL, ROR.UCL, p.value, p_adj, sig_label, direction) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print(n = Inf)

write_csv(simple_res, file.path(out_dir, "tier2_hlae_interaction_simple_ORs.csv"))
write_csv(ror_res, file.path(out_dir, "tier2_hlae_interaction_ROR_results.csv"))

message("\nDone. Outputs saved to: ", out_dir)
