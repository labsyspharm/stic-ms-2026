library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

# ── 1. Load & filter ─────────────────────────────────────────────────────────
df <- read_parquet(
  "C:/Users/tak415/OneDrive - Harvard University/Documents/CyCIF_Geomx_analysis_LSP_HMS/CYCIF_ORION_atlas_2.0/ORION-FISH/STIC_analysis/ROI_FISH_scoring_sc_with_FISH_group.parquet"
)

filt <- df %>%
  filter(
    `E-cadherinp` == TRUE,
    Incidental    == "Y",
    ROIcat        == "STIC"
  ) %>%
  mutate(
    both_LDH_Ki67 = (FISH_group == "Both") & (LDHp == TRUE) & (Ki67p == TRUE),
    both_LDH_only = (FISH_group == "Both") & (LDHp == TRUE) & (Ki67p == FALSE)
  )

# ── 2. One row per slide, two % columns ──────────────────────────────────────
freq_slide <- filt %>%
  group_by(slideName) %>%
  summarise(
    pct_LDH_Ki67 = mean(both_LDH_Ki67) * 100,
    pct_LDH_only = mean(both_LDH_only) * 100,
    .groups = "drop"
  )

# ── 3. Long format ───────────────────────────────────────────────────────────
freq_long <- freq_slide %>%
  pivot_longer(-slideName, names_to = "group", values_to = "pct") %>%
  mutate(
    group = recode(group,
      pct_LDH_Ki67 = "Both aneuploidy\nLDH+ Ki67+",
      pct_LDH_only = "Both aneuploidy\nLDH+ only"
    ),
    group = factor(group, levels = c(
      "Both aneuploidy\nLDH+ Ki67+",
      "Both aneuploidy\nLDH+ only"
    ))
  )

pal <- c(
  "Both aneuploidy\nLDH+ Ki67+" = "#E64B35",
  "Both aneuploidy\nLDH+ only"  = "#F4A261"
)

# ── 4. Plot ──────────────────────────────────────────────────────────────────
p <- ggplot(freq_long, aes(x = group, y = pct, colour = group)) +
  geom_boxplot(
    aes(fill = group),
    width = 0.35, outlier.shape = NA, alpha = 0.15, colour = "grey40"
  ) +
  geom_beeswarm(
    size = 3, cex = 2.5, shape = 21,
    aes(fill = group), colour = "white", stroke = 0.5
  ) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  stat_compare_means(
    method       = "wilcox.test",
    paired       = TRUE,
    comparisons  = list(c("Both aneuploidy\nLDH+ Ki67+", "Both aneuploidy\nLDH+ only")),
    label        = "p.format",
    tip.length   = 0.02,
    size         = 4
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(
    title    = "Both-aneuploid cells co-expressing LDH / Ki67 across slides",
    subtitle = "Filter: E-cadherin+ | Incidental Y | ROIcat == STIC\nEach dot = one slide  |  Denominator: all E-cad+ STIC cells per slide",
    x        = NULL,
    y        = "% of E-cad+ STIC cells per slide"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9),
    axis.text.x   = element_text(size = 11)
  )

ggsave("FISH_group_LDH_Ki67_beeswarm.pdf", p, width = 6, height = 5.5)
ggsave("FISH_group_LDH_Ki67_beeswarm.png", p, width = 6, height = 5.5, dpi = 150)
message("Saved.")
