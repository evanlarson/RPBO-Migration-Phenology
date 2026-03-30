options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

fig_dir <- "figures_jul_oct"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

is_valid_metric_for_phase <- function(phase_group, metric_label) {
  (phase_group == "departure" & metric_label != "first_detection") |
    (phase_group == "fall_arrival" & metric_label != "last_detection")
}

df <- read.csv("model_reduced_hypothesis_results.csv", stringsAsFactors = FALSE) %>%
  filter(
    is_valid_metric_for_phase(phase_group, metric_label),
    !is.na(temp_p_value_glm),
    !is.na(gam_temp_p_value),
    !is.na(delta_aic_gam_minus_glm),
    !is.na(gam_temp_edf)
  )

df <- df %>%
  mutate(
    sig_glm = !is.na(temp_p_fdr_global_glm) & temp_p_fdr_global_glm < 0.10,
    sig_gam = !is.na(temp_p_fdr_global_gam) & temp_p_fdr_global_gam < 0.10,
    sig_class = case_when(
      sig_glm & sig_gam ~ "Both GLM+GAM",
      sig_glm & !sig_gam ~ "GLM only",
      !sig_glm & sig_gam ~ "GAM only",
      TRUE ~ "Neither"
    ),
    sig_class = factor(sig_class, levels = c("Both GLM+GAM", "GLM only", "GAM only", "Neither")),
    neglog_glm = -log10(pmax(temp_p_value_glm, 1e-12)),
    neglog_gam = -log10(pmax(gam_temp_p_value, 1e-12)),
    gam_better = delta_aic_gam_minus_glm < -2
  )

sig_counts <- df %>%
  count(sig_class) %>%
  mutate(prop = n / sum(n))

count_lookup <- setNames(sig_counts$n, as.character(sig_counts$sig_class))
n_both <- ifelse("Both GLM+GAM" %in% names(count_lookup), count_lookup[["Both GLM+GAM"]], 0)
n_glm_only <- ifelse("GLM only" %in% names(count_lookup), count_lookup[["GLM only"]], 0)
n_gam_only <- ifelse("GAM only" %in% names(count_lookup), count_lookup[["GAM only"]], 0)

theme_set(theme_minimal(base_size = 11))
class_colors <- c(
  "Both GLM+GAM" = "#1b9e77",
  "GLM only" = "#d95f02",
  "GAM only" = "#7570b3",
  "Neither" = "#9aa0a6"
)

# Panel A: p-value agreement GLM vs GAM.
pA <- ggplot(df, aes(x = neglog_glm, y = neglog_gam, color = sig_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray50") +
  geom_point(size = 2.0, alpha = 0.82) +
  scale_color_manual(values = class_colors) +
  labs(
    title = "A. GLM vs GAM Significance Strength",
    x = expression(-log[10](p[GLM])),
    y = expression(-log[10](p[GAM])),
    color = ""
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# Panel B: AIC comparison.
pB <- ggplot(df, aes(x = delta_aic_gam_minus_glm, fill = gam_better)) +
  geom_histogram(binwidth = 1.5, color = "white", alpha = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55) +
  geom_vline(xintercept = -2, linetype = "dotted", linewidth = 0.55, color = "#7570b3") +
  scale_fill_manual(values = c(`TRUE` = "#7570b3", `FALSE` = "#bdbdbd"), labels = c(`TRUE` = "GAM better (ΔAIC<=-2)", `FALSE` = "No clear GAM gain")) +
  labs(
    title = "B. Model Fit Gain from GAM",
    x = "Delta AIC (GAM - GLM)",
    y = "Number of models",
    fill = ""
  ) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Panel C: GAM edf (nonlinearity strength).
pC <- ggplot(df, aes(x = gam_temp_edf)) +
  geom_histogram(binwidth = 0.2, fill = "#2a9d8f", color = "white", alpha = 0.9) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, color = "gray35") +
  labs(
    title = "C. GAM Linearity",
    x = "GAM edf for s(temp)",
    y = "Number of models"
  ) +
  theme(panel.grid.minor = element_blank())

# Panel D: FDR significance overlap.
pD <- ggplot(sig_counts, aes(x = sig_class, y = n, fill = sig_class)) +
  geom_col(width = 0.7, alpha = 0.92) +
  geom_text(aes(label = paste0(n, " (", round(prop * 100), "%)")), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "D. Significant Effects Detection by Method",
    x = "",
    y = "Number of models"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    panel.grid.minor = element_blank()
  )

out_file <- file.path(fig_dir, "figure_glm_vs_gam_comparison.png")
png(out_file, width = 4200, height = 3600, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

cat("Wrote:", out_file, "\n")
