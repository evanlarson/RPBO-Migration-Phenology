options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

fig_dir <- "figures_jul_oct"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

pretty_region <- function(x) {
  dplyr::recode(
    x,
    southern_coastal_bc = "Southern Coastal BC",
    central_coast = "Central Coast",
    south_interior = "South Interior",
    central_interior_sub_boreal = "Central Interior / Sub-Boreal",
    .default = x
  )
}

pretty_month <- function(x) {
  dplyr::recode(
    tolower(x),
    july = "July",
    august = "August",
    september = "September",
    october = "October",
    .default = x
  )
}

pretty_metric <- function(x) {
  dplyr::recode(
    x,
    first_detection = "First Detection",
    last_detection = "Last Detection",
    p10 = "P10",
    p50 = "P50",
    p90 = "P90",
    .default = x
  )
}

metrics_lookup <- data.frame(
  metric_code = c("first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"),
  metric_label = c("first_detection", "p10", "p50", "p90", "last_detection"),
  stringsAsFactors = FALSE
)

is_valid_metric_for_phase <- function(phase_group, metric_label) {
  (phase_group == "departure" & metric_label != "first_detection") |
    (phase_group == "fall_arrival" & metric_label != "last_detection")
}

month_var_lookup <- c(
  july = "july_mean_min_temp",
  august = "august_mean_min_temp",
  september = "september_mean_min_temp",
  october = "october_mean_min_temp"
)

safe_lm <- function(formula, data) {
  fit <- try(lm(formula, data = data), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

theme_set(theme_minimal(base_size = 11))

# ---------------------------------------------------------------------------
# Panel A: Regional warming signals (slopes with 95% CI).
# ---------------------------------------------------------------------------
weather_monthly <- read.csv("bc_weather_region_monthly_jul_oct.csv", stringsAsFactors = FALSE) %>%
  mutate(month_name = tolower(month_name))

weather_keys <- weather_monthly %>%
  distinct(region, month_name) %>%
  arrange(region, month_name)

weather_rows <- vector("list", nrow(weather_keys))
for (i in seq_len(nrow(weather_keys))) {
  rg <- weather_keys$region[i]
  mn <- weather_keys$month_name[i]
  dd <- weather_monthly %>%
    filter(region == rg, month_name == mn, !is.na(region_mean_min_temp), !is.na(year))

  out <- data.frame(
    region = rg,
    month_name = mn,
    n_years = nrow(dd),
    slope = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 8 && sd(dd$region_mean_min_temp) > 0) {
    fit <- safe_lm(region_mean_min_temp ~ year, dd)
    if (!is.null(fit)) {
      sm <- summary(fit)
      ci <- suppressWarnings(confint(fit, "year"))
      out$slope <- unname(coef(fit)["year"])
      out$conf_low <- ci[1]
      out$conf_high <- ci[2]
      out$p_value <- sm$coefficients["year", "Pr(>|t|)"]
    }
  }
  weather_rows[[i]] <- out
}

weather_trends <- bind_rows(weather_rows) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    sig = ifelse(!is.na(p_fdr) & p_fdr < 0.10, "FDR < 0.10", "Not significant"),
    region_label = pretty_region(region),
    month_label = pretty_month(month_name),
    label = paste(region_label, month_label, sep = " | ")
  ) %>%
  arrange(slope)

pA <- ggplot(weather_trends, aes(x = slope, y = reorder(label, slope), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray45") +
  geom_errorbar(
    aes(y = reorder(label, slope), xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.7
  ) +
  geom_point(size = 2.1) +
  scale_color_manual(values = c("FDR < 0.10" = "#d95f02", "Not significant" = "#6c757d")) +
  labs(
    title = "A. Regional Temperature Trends During Migration Months (July-October)",
    x = "Temperature Trend (deg C per year)",
    y = "",
    color = ""
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ---------------------------------------------------------------------------
# Panel B: Species timing trends over years (no broad species-wide shift).
# ---------------------------------------------------------------------------
phen <- read.csv("rpbo_focal_species_phenology_jul_oct.csv", stringsAsFactors = FALSE)

phen_long <- phen %>%
  select(
    year, SPECIES, SPNAME, phase_group,
    first_julian_day, q10_julian_day, q50_julian_day, q90_julian_day, last_julian_day
  ) %>%
  pivot_longer(
    cols = c(first_julian_day, q10_julian_day, q50_julian_day, q90_julian_day, last_julian_day),
    names_to = "metric_code",
    values_to = "timing_julian_day"
  ) %>%
  left_join(metrics_lookup, by = "metric_code")

phen_keys <- phen_long %>%
  distinct(SPECIES, SPNAME, phase_group, metric_label) %>%
  arrange(phase_group, SPNAME, metric_label)

phen_rows <- vector("list", nrow(phen_keys))
for (i in seq_len(nrow(phen_keys))) {
  sp <- phen_keys$SPECIES[i]
  mt <- phen_keys$metric_label[i]
  dd <- phen_long %>%
    filter(SPECIES == sp, metric_label == mt, !is.na(timing_julian_day), !is.na(year))

  out <- data.frame(
    SPECIES = sp,
    SPNAME = phen_keys$SPNAME[i],
    phase_group = phen_keys$phase_group[i],
    metric_label = mt,
    n_years = nrow(dd),
    slope = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 8 && sd(dd$timing_julian_day) > 0) {
    fit <- safe_lm(timing_julian_day ~ year, dd)
    if (!is.null(fit)) {
      sm <- summary(fit)
      ci <- suppressWarnings(confint(fit, "year"))
      out$slope <- unname(coef(fit)["year"])
      out$conf_low <- ci[1]
      out$conf_high <- ci[2]
      out$p_value <- sm$coefficients["year", "Pr(>|t|)"]
    }
  }
  phen_rows[[i]] <- out
}

phen_trends <- bind_rows(phen_rows) %>%
  filter(is_valid_metric_for_phase(phase_group, metric_label)) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    sig = ifelse(!is.na(p_fdr) & p_fdr < 0.10, "FDR < 0.10", "Not significant"),
    metric_short = dplyr::recode(
      metric_label,
      first_detection = "first",
      last_detection = "last",
      .default = metric_label
    ),
    label = paste(SPNAME, metric_short, sep = " | ")
  ) %>%
  arrange(slope)

pB <- ggplot(phen_trends, aes(x = slope, y = reorder(label, slope), color = phase_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray45") +
  geom_errorbar(
    aes(y = reorder(label, slope), xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(size = 1.9) +
  scale_color_manual(
    values = c(departure = "#1b9e77", fall_arrival = "#d95f02"),
    labels = c(departure = "Departure", fall_arrival = "Arrival")
  ) +
  labs(
    title = "B. Species-Specific Phenology Trends by Year",
    x = "Timing Trend (days per year)",
    y = "",
    color = "Group"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# Panel C: Most robust temperature-timing relationships from reduced tests.
# ---------------------------------------------------------------------------
reduced <- read.csv("model_reduced_hypothesis_results.csv", stringsAsFactors = FALSE)
weather_wide <- read.csv("bc_weather_region_monthly_wide_jul_oct.csv", stringsAsFactors = FALSE)

robust <- reduced %>%
  filter(
    temp_p_fdr_global_glm < 0.10 |
      temp_p_fdr_global_permutation < 0.10 |
      temp_p_fdr_global_gam < 0.10
  ) %>%
  filter(is_valid_metric_for_phase(phase_group, metric_label)) %>%
  arrange(temp_p_fdr_global_glm, temp_p_fdr_global_permutation, temp_p_fdr_global_gam, temp_p_value_glm) %>%
  select(SPECIES, SPNAME, phase_group, metric_label, region, temp_month, temp_slope_days_per_degC, temp_p_fdr_global_glm, temp_p_fdr_global_permutation, temp_p_fdr_global_gam) %>%
  distinct()

build_relationship_data <- function(r) {
  d_sp <- phen_long %>%
    filter(SPECIES == r$SPECIES, metric_label == r$metric_label) %>%
    select(year, timing_julian_day)
  d_rg <- weather_wide %>%
    filter(region == r$region)
  mvar <- month_var_lookup[[r$temp_month]]
  if (is.null(mvar) || !(mvar %in% names(d_rg))) return(NULL)

  dd <- inner_join(d_sp, d_rg %>% select(year, all_of(mvar)), by = "year") %>%
    rename(temp_value = all_of(mvar)) %>%
    filter(!is.na(timing_julian_day), !is.na(temp_value))

  if (!nrow(dd)) return(NULL)

  dd$panel <- sprintf(
    "%s | %s | %s temp\n%s (%+.2f d/deg C)",
    r$SPNAME,
    pretty_region(r$region),
    pretty_month(r$temp_month),
    pretty_metric(r$metric_label),
    r$temp_slope_days_per_degC
  )
  dd
}

rel_list <- lapply(seq_len(nrow(robust)), function(i) build_relationship_data(robust[i, ]))
rel_data <- bind_rows(rel_list)

pC <- ggplot(rel_data, aes(x = temp_value, y = timing_julian_day)) +
  geom_point(size = 1.8, alpha = 0.85, color = "#2d3748") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "#e63946") +
  facet_wrap(~panel, scales = "free", ncol = 2) +
  labs(
    title = "C. Strongest Temperature-Timing Associations",
    x = "Regional Monthly Mean Minimum Temperature (deg C)",
    y = "Timing Metric (Julian day)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8.8)
  )

# ---------------------------------------------------------------------------
# Compose into a single concept figure (3 stacked panels).
# ---------------------------------------------------------------------------
out_file <- file.path(fig_dir, "figure_final_concept_summary.png")
png(out_file, width = 4200, height = 5400, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 3, ncol = 1,
  heights = unit(c(1.0, 1.2, 1.4), "null")
)))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pC, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()

cat("Wrote:", out_file, "\n")
