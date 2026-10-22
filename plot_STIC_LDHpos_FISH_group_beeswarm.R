library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# ── 1. Load ───────────────────────────────────────────────────────────────────
df <- read_parquet(
  "C:/Users/tak415/OneDrive - Harvard University/Documents/CyCIF_Geomx_analysis_LSP_HMS/CYCIF_ORION_atlas_2.0/ORION-FISH/STIC_analysis/ROI_FISH_scoring_sc_with_FISH_group.parquet"
)

# ── 2. Filter: STIC, incidental, E-cad+, LDH+ ────────────────────────────────
stic_ldh <- df %>%
  filter(
    `E-cadherinp` == TRUE,
    Incidental    == "Y",
    ROIcat        == "STIC",
    LDHp          == TRUE
  )

# ── 3. Per slide: % of LDH+ cells in each FISH_group ─────────────────────────
# Denominator = total LDH+ cells per slide
stic_freq <- stic_ldh %>%
  group_by(slideName) %>%
  summarise(
    n_ldh       = n(),
    pct_Both    = mean(FISH_group == "Both")       * 100,
    pct_MYC     = mean(FISH_group == "MYC only")   * 100,
    pct_CCNE1   = mean(FISH_group == "CCNE1 only") * 100,
    pct_Diploid = mean(FISH_group == "Diploid")    * 100,
    .groups     = "drop"
  )

cat("=== STIC-only per-slide summary ===\n")
print(stic_freq)

# ── 4. Long format for plotting ──────────────────────────────────────────────
stic_long <- stic_freq %>%
  pivot_longer(
    cols      = starts_with("pct_"),
    names_to  = "FISH_group",
    values_to = "pct"
  ) %>%
  mutate(
    FISH_group = recode(FISH_group,
      pct_Both    = "Both aneuploid",
      pct_MYC     = "MYC only",
      pct_CCNE1   = "CCNE1 only",
      pct_Diploid = "Diploid"
    ),
    FISH_group = factor(FISH_group,
      levels = c("Both aneuploid", "MYC only", "CCNE1 only", "Diploid"))
  )

fish_pal4 <- c(
  "Both aneuploid" = "#E64B35",
  "MYC only"       = "#F4A261",
  "CCNE1 only"     = "#4DBBD5",
  "Diploid"        = "#B0B0B0"
)

# ── 5. Beeswarm: % of LDH+ cells per slide, by FISH group (STIC only) ───────
p <- ggplot(stic_long, aes(x = FISH_group, y = pct, colour = FISH_group, fill = FISH_group)) +
  geom_boxplot(
    width = 0.4, outlier.shape = NA, alpha = 0.15, colour = "grey40"
  ) +
  geom_beeswarm(
    size = 3, cex = 2.5, shape = 21, colour = "white", stroke = 0.5
  ) +
  scale_colour_manual(values = fish_pal4, guide = "none") +
  scale_fill_manual(values = fish_pal4, guide = "none") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title    = "FISH_group composition of LDH+ cells (STIC only)",
    subtitle = "Filter: E-cad+ | Incidental Y | ROIcat == STIC\nEach dot = one slide  |  Denominator: LDH+ cells per slide",
    x        = NULL,
    y        = "% of LDH+ cells per slide"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9),
    axis.text.x   = element_text(size = 11)
  )

ggsave("STIC_LDHpos_FISH_group_beeswarm.pdf", p, width = 7, height = 5.5)
ggsave("STIC_LDHpos_FISH_group_beeswarm.png", p, width = 7, height = 5.5, dpi = 150)
message("Saved: STIC_LDHpos_FISH_group_beeswarm.pdf / .png")
