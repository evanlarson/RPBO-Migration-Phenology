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

month_levels <- c("july", "august", "september", "october")
region_levels <- c(
  "southern_coastal_bc",
  "central_coast",
  "south_interior",
  "central_interior_sub_boreal"
)

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
    x,
    july = "July",
    august = "August",
    september = "September",
    october = "October",
    .default = x
  )
}

theme_set(theme_minimal(base_size = 12))

# Shared equal-interval date axis (every 14 days across the banding season).
equal_tick_dates <- seq(as.Date("2001-07-21"), as.Date("2001-10-18"), by = "14 days")
equal_tick_days <- as.integer(format(equal_tick_dates, "%j"))
equal_tick_labels <- format(equal_tick_dates, "%b %d")

# Input tables.
reg_month <- read.csv("bc_weather_region_monthly_jul_oct.csv", stringsAsFactors = FALSE) %>%
  mutate(
    region = factor(region, levels = region_levels, labels = pretty_region(region_levels)),
    month_name = factor(month_name, levels = month_levels, labels = pretty_month(month_levels))
  )

st_month <- read.csv("bc_weather_station_monthly_jul_oct.csv", stringsAsFactors = FALSE) %>%
  mutate(
    region_label = pretty_region(region),
    station_panel = paste0(station, "\n(", region_label, ")"),
    month_name = factor(month_name, levels = month_levels, labels = pretty_month(month_levels)),
    station_panel = factor(station_panel, levels = unique(station_panel[order(region_label, station)]))
  )

# 1) Composite station monthly temperatures:
# 12 station subsections; each has Jul/Aug/Sep/Oct scatter + trendline.
p_station <- ggplot(st_month, aes(x = year, y = mean_min_temp, color = month_name)) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~station_panel, ncol = 4, scales = "free_y") +
  scale_color_brewer(palette = "Dark2", name = "Month") +
  labs(
    title = "Station Monthly Minimum Temperatures (July-October)",
    subtitle = "Each panel is one station. Points are monthly values; lines are month-specific linear trends.",
    x = "Year",
    y = "Mean Minimum Temperature (deg C)",
    color = "Month"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(
  filename = file.path(fig_dir, "figure_station_monthly_composite_scatter_trends.png"),
  plot = p_station,
  width = 16,
  height = 11,
  dpi = 300
)

# 2) Regional monthly temperatures:
# 4 region panels with month-colored scatter + trendline.
region_trend_stats <- reg_month %>%
  group_by(region, month_name) %>%
  summarize(
    slope = unname(coef(lm(region_mean_min_temp ~ year))[["year"]]),
    intercept = unname(coef(lm(region_mean_min_temp ~ year))[["(Intercept)"]]),
    p_value = summary(lm(region_mean_min_temp ~ year))$coefficients["year", "Pr(>|t|)"],
    x_start = min(year, na.rm = TRUE),
    x_end = max(year, na.rm = TRUE),
    .groups = "drop"
  )

region_sig_markers <- region_trend_stats %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  mutate(
    x_star = x_start + 0.45,
    y_line = intercept + slope * x_star,
    sig_label = "*"
  ) %>%
  group_by(region) %>%
  arrange(desc(y_line), .by_group = TRUE) %>%
  group_modify(~{
    d <- .x
    if (nrow(d) == 0) return(d)
    min_gap <- 0.35
    y_adj <- d$y_line
    if (length(y_adj) > 1) {
      for (i in 2:length(y_adj)) {
        if ((y_adj[i - 1] - y_adj[i]) < min_gap) {
          y_adj[i] <- y_adj[i - 1] - min_gap
        }
      }
    }
    d$y_star <- y_adj
    d
  }) %>%
  ungroup()

p_region <- ggplot(reg_month, aes(x = year, y = region_mean_min_temp, color = month_name)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.0) +
  geom_text(
    data = region_sig_markers,
    aes(x = x_star, y = y_star, label = sig_label),
    inherit.aes = FALSE,
    color = "black",
    size = 8.4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~region, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c(
      "July" = "#f1c40f",
      "August" = "#e67e22",
      "September" = "#2ca25f",
      "October" = "#2b8cbe"
    ),
    name = "Month"
  ) +
  labs(
    title = "Regional Monthly Minimum Temperatures (July-October)",
    x = "Year",
    y = "Mean Minimum Temperature (deg C)"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_regional_monthly_composite_scatter_trends.png"),
  plot = p_region,
  width = 13,
  height = 9,
  dpi = 300
)

# 3) Median passage date by species over years.
phen <- read.csv("rpbo_focal_species_phenology_jul_oct.csv", stringsAsFactors = FALSE) %>%
  mutate(
    phase_group = factor(phase_group, levels = c("departure", "fall_arrival")),
    SPNAME = factor(SPNAME, levels = unique(SPNAME[order(phase_group, SPNAME)]))
  )

p_median <- ggplot(phen, aes(x = year, y = q50_julian_day, color = phase_group)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~SPNAME, ncol = 3) +
  scale_color_manual(
    values = c(departure = "#1b9e77", fall_arrival = "#d95f02"),
    labels = c(departure = "Departure species", fall_arrival = "Fall-arrival species")
  ) +
  labs(
    title = "Median Passage Date (q50) by Species Over Years",
    subtitle = "July-October filtered data, 2005-2024",
    x = "Year",
    y = "Median Passage Date (Julian Day)",
    color = "Species Group"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_species_median_passage_vs_year.png"),
  plot = p_median,
  width = 12,
  height = 8,
  dpi = 300
)

# 4) Multi-metric timing plot for fall-arrival species (first, p10, p50, p90).
arrive_metrics <- phen %>%
  filter(phase_group == "fall_arrival") %>%
  select(year, SPECIES, SPNAME, first_julian_day, q10_julian_day, q50_julian_day, q90_julian_day) %>%
  pivot_longer(
    cols = c(first_julian_day, q10_julian_day, q50_julian_day, q90_julian_day),
    names_to = "metric",
    values_to = "metric_day"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day"),
      labels = c("First detection", "p10", "p50", "p90")
    )
  )

p_arrive_multi <- ggplot(arrive_metrics, aes(x = year, y = metric_day, color = metric)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~SPNAME, ncol = 3) +
  scale_color_manual(
    values = c(
      "First detection" = "#e41a1c",
      "p10" = "#377eb8",
      "p50" = "#4daf4a",
      "p90" = "#984ea3"
    ),
    name = "Metric"
  ) +
  scale_y_continuous(
    breaks = equal_tick_days,
    labels = equal_tick_labels
  ) +
  labs(
    title = "Phenology Metrics for Arrival Species",
    x = "Year",
    y = "Calendar Date"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_arrival_multimetric_dates_by_year.png"),
  plot = p_arrive_multi,
  width = 12,
  height = 8,
  dpi = 300
)

# 5) Multi-metric timing plot for departure species (p10, p50, p90, last).
depart_metrics <- phen %>%
  filter(phase_group == "departure") %>%
  select(year, SPECIES, SPNAME, q10_julian_day, q50_julian_day, q90_julian_day, last_julian_day) %>%
  pivot_longer(
    cols = c(q10_julian_day, q50_julian_day, q90_julian_day, last_julian_day),
    names_to = "metric",
    values_to = "metric_day"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"),
      labels = c("p10", "p50", "p90", "Last detection")
    )
  )

p_depart_multi <- ggplot(depart_metrics, aes(x = year, y = metric_day, color = metric)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~SPNAME, ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Metric") +
  scale_y_continuous(
    breaks = equal_tick_days,
    labels = equal_tick_labels
  ) +
  labs(
    title = "Phenology Metrics for Departure Species",
    x = "Year",
    y = "Calendar Date"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_departure_multimetric_dates_by_year.png"),
  plot = p_depart_multi,
  width = 12,
  height = 8,
  dpi = 300
)

# 5b) Combined arrival + departure phenology metric figure.
combined_multimetric_file <- file.path(fig_dir, "figure_arrival_departure_multimetric_dates_by_year.png")
png(combined_multimetric_file, width = 3600, height = 5000, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
print(
  p_depart_multi + theme(plot.margin = margin(6, 8, 3, 8)),
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
)
print(
  p_arrive_multi + theme(plot.margin = margin(3, 8, 6, 8)),
  vp = viewport(layout.pos.row = 2, layout.pos.col = 1)
)
dev.off()

# 6) First capture (arrival species) and last capture (departure species) by year.
first_last <- read.csv("rpbo_focal_species_phenology_jul_oct.csv", stringsAsFactors = FALSE) %>%
  mutate(
    capture_metric_day = ifelse(phase_group == "fall_arrival", first_julian_day, last_julian_day),
    capture_metric_label = ifelse(
      phase_group == "fall_arrival",
      "First capture in season",
      "Last capture in season"
    ),
    phase_group = factor(phase_group, levels = c("departure", "fall_arrival")),
    SPNAME = factor(SPNAME, levels = unique(SPNAME[order(phase_group, SPNAME)]))
  )

p_first_last <- ggplot(first_last, aes(x = year, y = capture_metric_day, color = phase_group)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~SPNAME, ncol = 3) +
  scale_color_manual(
    values = c(departure = "#1b9e77", fall_arrival = "#d95f02"),
    labels = c(
      departure = "Departure species (last capture)",
      fall_arrival = "Fall-arrival species (first capture)"
    )
  ) +
  labs(
    title = "First/Last Capture Timing by Year",
    subtitle = "Arrival species: first capture date each year. Departure species: last capture date each year.",
    x = "Year",
    y = "Capture Timing (Julian Day)",
    color = "Metric"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_first_last_capture_by_year.png"),
  plot = p_first_last,
  width = 12,
  height = 8,
  dpi = 300
)

# 7) Seasonal banding profile per species (rise/fall through season).
species_groups <- read.csv("focal_species_groups_jul_oct.csv", stringsAsFactors = FALSE) %>%
  select(SPECIES, SPNAME, phase_group)

rpbo_all <- read.csv("rpbo_clean_standardized.csv", stringsAsFactors = FALSE) %>%
  mutate(
    date = as.Date(date_iso),
    month = as.integer(format(date, "%m"))
  ) %>%
  filter(!is.na(month), month %in% c(7, 8, 9, 10), !is.na(year_std), !is.na(julian_day), !is.na(Banded)) %>%
  inner_join(species_groups, by = c("SPECIES", "SPNAME"))

daily_species_year <- rpbo_all %>%
  group_by(SPECIES, SPNAME, phase_group, year_std, julian_day) %>%
  summarize(banded_day = sum(Banded, na.rm = TRUE), .groups = "drop")

seasonal_profile <- daily_species_year %>%
  group_by(SPECIES, SPNAME, phase_group, julian_day) %>%
  summarize(
    mean_daily_banded = mean(banded_day, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(phase_group = factor(phase_group, levels = c("departure", "fall_arrival")))

seasonal_profile_plot <- seasonal_profile %>%
  filter(julian_day >= 202, julian_day <= 291)

phase_labels <- c(
  departure = "Departures",
  fall_arrival = "Arrivals"
)

# Explicit warm/cool palette (3 species each) with high contrast.
warm_cols <- c(
  "#D62828", # red
  "#F28E2B", # orange
  "#C8A100"  # yellow (darker for visibility)
)
cool_cols <- c(
  "#1B9E77", # green
  "#56B4E9", # light blue
  "#1F3A93"  # dark blue
)

dep_names <- seasonal_profile %>%
  filter(phase_group == "departure") %>%
  distinct(SPNAME) %>%
  arrange(SPNAME) %>%
  pull(SPNAME)

arr_names <- seasonal_profile %>%
  filter(phase_group == "fall_arrival") %>%
  distinct(SPNAME) %>%
  arrange(SPNAME) %>%
  pull(SPNAME)

species_colors <- c(
  setNames(warm_cols[seq_along(dep_names)], dep_names),
  setNames(cool_cols[seq_along(arr_names)], arr_names)
)
legend_breaks <- c(dep_names, arr_names)

p_seasonal <- ggplot(seasonal_profile_plot, aes(x = julian_day, y = mean_daily_banded, color = SPNAME)) +
  geom_point(size = 1.0, alpha = 0.45) +
  geom_smooth(method = "loess", se = FALSE, span = 0.25, linewidth = 1.1) +
  facet_wrap(
    ~phase_group,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(phase_labels)
  ) +
  scale_x_continuous(
    breaks = equal_tick_days,
    labels = equal_tick_labels,
    limits = c(202, 291)
  ) +
  scale_color_manual(values = species_colors, breaks = legend_breaks) +
  labs(
    title = "Migration Phenology at Rocky Point Bird Observatory",
    x = "Calendar Date",
    y = "Mean Banded per Day",
    color = "Species"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "figure_species_seasonal_profile.png"),
  plot = p_seasonal,
  width = 12,
  height = 9,
  dpi = 300
)

# 8) Moving-average peaks (3-, 5-, 7-day) for each species.
all_days <- data.frame(julian_day = 182:304)
species_levels <- unique(seasonal_profile[, c("SPECIES", "SPNAME", "phase_group")])

moving_rows <- list()
peak_rows <- list()

for (i in seq_len(nrow(species_levels))) {
  sp <- species_levels$SPECIES[i]
  sp_name <- species_levels$SPNAME[i]
  sp_phase <- species_levels$phase_group[i]

  d <- seasonal_profile %>%
    filter(SPECIES == sp) %>%
    select(julian_day, mean_daily_banded)

  d_full <- all_days %>%
    left_join(d, by = "julian_day") %>%
    mutate(mean_daily_banded = ifelse(is.na(mean_daily_banded), 0, mean_daily_banded))

  for (k in c(3, 5, 7)) {
    ma <- as.numeric(stats::filter(d_full$mean_daily_banded, rep(1 / k, k), sides = 2))
    tmp <- data.frame(
      SPECIES = sp,
      SPNAME = sp_name,
      phase_group = sp_phase,
      julian_day = d_full$julian_day,
      window = paste0(k, "-day"),
      moving_avg = ma,
      stringsAsFactors = FALSE
    )
    moving_rows[[length(moving_rows) + 1]] <- tmp

    if (all(is.na(ma))) next
    peak_idx <- which(ma == max(ma, na.rm = TRUE))[1]
    peak_rows[[length(peak_rows) + 1]] <- data.frame(
      SPECIES = sp,
      SPNAME = sp_name,
      phase_group = sp_phase,
      window = paste0(k, "-day"),
      peak_julian_day = d_full$julian_day[peak_idx],
      peak_value = ma[peak_idx],
      stringsAsFactors = FALSE
    )
  }
}

moving_df <- bind_rows(moving_rows) %>%
  mutate(
    window = factor(window, levels = c("3-day", "5-day", "7-day")),
    SPNAME = factor(SPNAME, levels = unique(SPNAME[order(phase_group, SPNAME)]))
  )

peak_df <- bind_rows(peak_rows) %>%
  mutate(
    window = factor(window, levels = c("3-day", "5-day", "7-day")),
    SPNAME = factor(SPNAME, levels = levels(moving_df$SPNAME))
  )

p_moving <- ggplot(moving_df, aes(x = julian_day, y = moving_avg, color = window)) +
  geom_line(linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
  geom_point(
    data = peak_df,
    aes(x = peak_julian_day, y = peak_value, fill = window),
    shape = 21,
    size = 2.3,
    stroke = 0.6,
    color = "black"
  ) +
  facet_wrap(~SPNAME, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("3-day" = "#1b9e77", "5-day" = "#d95f02", "7-day" = "#7570b3")) +
  scale_fill_manual(values = c("3-day" = "#1b9e77", "5-day" = "#d95f02", "7-day" = "#7570b3")) +
  labs(
    title = "Moving-Average Migration Peaks by Species",
    subtitle = "Lines are 3-, 5-, and 7-day moving averages of mean daily banding; points mark each line's peak day.",
    x = "Julian Day",
    y = "Moving-Average Banded per Day",
    color = "Window",
    fill = "Window"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "figure_species_moving_average_peaks.png"),
  plot = p_moving,
  width = 13,
  height = 9,
  dpi = 300
)

cat("Saved figures to:", fig_dir, "\n")
for (f in list.files(fig_dir, pattern = "^figure_(station|regional|species_median|arrival_multimetric|departure_multimetric|first_last|species_seasonal|species_moving)", full.names = FALSE)) {
  cat(" -", f, "\n")
}
