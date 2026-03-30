options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(mgcv)
})

# Inputs
phenology_file <- "rpbo_focal_species_phenology_jul_oct.csv"
weather_monthly_file <- "bc_weather_region_monthly_jul_oct.csv"
weather_wide_file <- "bc_weather_region_monthly_wide_jul_oct.csv"
net_hours_file <- "nethours.xlsx"

# Outputs
out_weather_trend <- "model_weather_warming_trends_glm.csv"
out_phenology_trend <- "model_phenology_year_trends_glm.csv"
out_glm_assoc <- "model_phenology_temperature_associations_glm.csv"
out_glm_summary <- "model_phenology_temperature_associations_glm_summary.csv"
out_glm_sig <- "model_phenology_temperature_associations_glm_significant.csv"
out_gam_assoc <- "model_phenology_temperature_associations_gam.csv"
out_gam_sig <- "model_phenology_temperature_associations_gam_significant.csv"

metrics_lookup <- data.frame(
  metric_code = c("first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"),
  metric_label = c("first_detection", "p10", "p50", "p90", "last_detection"),
  stringsAsFactors = FALSE
)

is_valid_metric_for_phase <- function(phase_group, metric_label) {
  (phase_group == "departure" & metric_label != "first_detection") |
    (phase_group == "fall_arrival" & metric_label != "last_detection")
}

month_vars <- c("july_mean_min_temp", "august_mean_min_temp", "september_mean_min_temp", "october_mean_min_temp")
month_labels <- c("july", "august", "september", "october")
names(month_labels) <- month_vars

safe_lm <- function(formula, data) {
  fit <- try(lm(formula, data = data), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

safe_glm <- function(formula, data) {
  fit <- try(glm(formula, data = data, family = gaussian()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

safe_gam <- function(formula, data) {
  fit <- try(gam(formula, data = data, method = "REML"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

load_net_hours <- function(path) {
  if (!file.exists(path)) {
    stop(paste("Net-hours file not found:", path))
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl is required to read nethours.xlsx")
    }
    d <- as.data.frame(readxl::read_excel(path))
  } else {
    d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  nms <- tolower(gsub("[^a-z0-9]+", "_", names(d)))
  names(d) <- nms

  year_col <- names(d)[nms == "year"][1]
  if (is.na(year_col)) {
    stop("Net-hours file must include a 'year' column.")
  }

  net_col_candidates <- names(d)[grepl("net", nms) & grepl("hour", nms)]
  if (!length(net_col_candidates)) {
    stop("Net-hours file must include a net-hours column (e.g., 'net hours').")
  }
  net_col <- net_col_candidates[1]

  out <- data.frame(
    year = suppressWarnings(as.integer(d[[year_col]])),
    net_hours = suppressWarnings(as.numeric(d[[net_col]])),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(year), !is.na(net_hours), net_hours > 0) %>%
    group_by(year) %>%
    summarize(net_hours = mean(net_hours, na.rm = TRUE), .groups = "drop") %>%
    mutate(net_hours_100 = net_hours / 100)

  out
}

interpret_direction <- function(slope) {
  ifelse(
    is.na(slope),
    NA_character_,
    ifelse(slope < 0, "warmer_or_later_covariate -> earlier_timing", "warmer_or_later_covariate -> later_timing")
  )
}

# Load data.
phen <- read.csv(phenology_file, stringsAsFactors = FALSE)
weather_monthly <- read.csv(weather_monthly_file, stringsAsFactors = FALSE)
weather_wide <- read.csv(weather_wide_file, stringsAsFactors = FALSE)
net_hours <- load_net_hours(net_hours_file)

required_phen <- c(
  "year", "SPECIES", "SPNAME", "phase_group",
  "first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"
)
missing_phen <- setdiff(required_phen, names(phen))
if (length(missing_phen) > 0) {
  stop(paste("Missing phenology columns:", paste(missing_phen, collapse = ", ")))
}

required_wm <- c("region", "year", "month_name", "region_mean_min_temp")
missing_wm <- setdiff(required_wm, names(weather_monthly))
if (length(missing_wm) > 0) {
  stop(paste("Missing weather monthly columns:", paste(missing_wm, collapse = ", ")))
}

required_ww <- c("region", "year", month_vars)
missing_ww <- setdiff(required_ww, names(weather_wide))
if (length(missing_ww) > 0) {
  stop(paste("Missing weather wide columns:", paste(missing_ww, collapse = ", ")))
}

# Long phenology table for first / p10 / p50 / p90 / last.
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

phen_long <- phen_long %>%
  filter(is_valid_metric_for_phase(phase_group, metric_label))

phen_long <- phen_long %>%
  left_join(net_hours, by = "year")

# 1) Regional temperature trend over years (GLM/LM): temp ~ year
weather_monthly <- weather_monthly %>%
  mutate(month_name = tolower(month_name)) %>%
  filter(month_name %in% month_labels)

wm_keys <- weather_monthly %>% distinct(region, month_name)
wm_rows <- vector("list", nrow(wm_keys))

for (i in seq_len(nrow(wm_keys))) {
  rg <- wm_keys$region[i]
  mn <- wm_keys$month_name[i]
  dd <- weather_monthly %>% filter(region == rg, month_name == mn, !is.na(region_mean_min_temp), !is.na(year))

  out <- data.frame(
    region = rg,
    month_name = mn,
    n_years = nrow(dd),
    slope_degC_per_year = NA_real_,
    p_value = NA_real_,
    r_squared = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 8 && stats::sd(dd$region_mean_min_temp) > 0) {
    fit <- safe_lm(region_mean_min_temp ~ year, dd)
    if (!is.null(fit)) {
      sm <- summary(fit)
      out$slope_degC_per_year <- unname(coef(fit)[2])
      out$p_value <- sm$coefficients[2, 4]
      out$r_squared <- sm$r.squared
    }
  }
  wm_rows[[i]] <- out
}

weather_trend <- bind_rows(wm_rows) %>%
  mutate(
    p_value_fdr_global = p.adjust(p_value, method = "BH"),
    trend_direction = ifelse(is.na(slope_degC_per_year), NA_character_, ifelse(slope_degC_per_year < 0, "cooling", "warming"))
  ) %>%
  arrange(p_value)

write.csv(weather_trend, out_weather_trend, row.names = FALSE, na = "")

# 2) Phenology trend over years, effort-adjusted: timing ~ year + net_hours_100
pm_keys <- phen_long %>% distinct(SPECIES, SPNAME, phase_group, metric_label)
pm_rows <- vector("list", nrow(pm_keys))

for (i in seq_len(nrow(pm_keys))) {
  sp <- pm_keys$SPECIES[i]
  mt <- pm_keys$metric_label[i]
  dd <- phen_long %>%
    filter(
      SPECIES == sp,
      metric_label == mt,
      !is.na(timing_julian_day),
      !is.na(year),
      !is.na(net_hours_100)
    )

  out <- data.frame(
    SPECIES = sp,
    SPNAME = pm_keys$SPNAME[i],
    phase_group = pm_keys$phase_group[i],
    metric_label = mt,
    n_years = nrow(dd),
    year_slope_days_per_year = NA_real_,
    year_p_value = NA_real_,
    net_hours_slope_days_per_100h = NA_real_,
    net_hours_p_value = NA_real_,
    r_squared = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 8 && stats::sd(dd$timing_julian_day) > 0 && stats::sd(dd$net_hours_100) > 0) {
    fit <- safe_lm(timing_julian_day ~ year + net_hours_100, dd)
    if (!is.null(fit)) {
      sm <- summary(fit)
      if ("year" %in% rownames(sm$coefficients)) {
        out$year_slope_days_per_year <- unname(coef(fit)["year"])
        out$year_p_value <- sm$coefficients["year", 4]
      }
      if ("net_hours_100" %in% rownames(sm$coefficients)) {
        out$net_hours_slope_days_per_100h <- unname(coef(fit)["net_hours_100"])
        out$net_hours_p_value <- sm$coefficients["net_hours_100", 4]
      }
      out$r_squared <- sm$r.squared
    }
  }
  pm_rows[[i]] <- out
}

phenology_trend <- bind_rows(pm_rows) %>%
  mutate(
    p_value = year_p_value,
    slope_days_per_year = year_slope_days_per_year,
    p_value_fdr_global = p.adjust(year_p_value, method = "BH"),
    trend_direction = ifelse(is.na(year_slope_days_per_year), NA_character_, ifelse(year_slope_days_per_year < 0, "earlier_over_time", "later_over_time"))
  ) %>%
  arrange(year_p_value)

write.csv(phenology_trend, out_phenology_trend, row.names = FALSE, na = "")

# 3) Main association model (GLM Gaussian), effort-adjusted: timing ~ year + net_hours_100 + monthly_temp
glm_rows <- list()
k <- 0
for (i in seq_len(nrow(pm_keys))) {
  sp <- pm_keys$SPECIES[i]
  spname <- pm_keys$SPNAME[i]
  phase <- pm_keys$phase_group[i]
  mt <- pm_keys$metric_label[i]

  d_sp <- phen_long %>%
    filter(SPECIES == sp, metric_label == mt) %>%
    select(year, timing_julian_day, net_hours_100)

  regions <- unique(weather_wide$region)
  for (rg in regions) {
    d_rg <- weather_wide %>% filter(region == rg)
    merged <- inner_join(d_sp, d_rg, by = "year")

    for (mv in month_vars) {
      temp_label <- month_labels[[mv]]
      dd <- merged %>%
        transmute(
          year = year,
          timing_julian_day = timing_julian_day,
          net_hours_100 = net_hours_100,
          temp = .data[[mv]]
        ) %>%
        filter(!is.na(timing_julian_day), !is.na(temp), !is.na(net_hours_100), !is.na(year))

      k <- k + 1
      out <- data.frame(
        SPECIES = sp,
        SPNAME = spname,
        phase_group = phase,
        metric_label = mt,
        region = rg,
        temp_month = temp_label,
        n_years = nrow(dd),
        temp_slope_days_per_degC = NA_real_,
        temp_p_value = NA_real_,
        net_hours_slope_days_per_100h = NA_real_,
        net_hours_p_value = NA_real_,
        year_slope_days_per_year = NA_real_,
        year_p_value = NA_real_,
        model_r_squared = NA_real_,
        model_aic = NA_real_,
        stringsAsFactors = FALSE
      )

      if (nrow(dd) >= 10 &&
          stats::sd(dd$timing_julian_day) > 0 &&
          stats::sd(dd$net_hours_100) > 0 &&
          stats::sd(dd$temp) > 0) {
        fit <- safe_glm(timing_julian_day ~ year + net_hours_100 + temp, dd)
        if (!is.null(fit)) {
          sm <- summary(fit)
          cf <- sm$coefficients
          if ("temp" %in% rownames(cf)) {
            out$temp_slope_days_per_degC <- unname(coef(fit)["temp"])
            out$temp_p_value <- cf["temp", 4]
          }
          if ("year" %in% rownames(cf)) {
            out$year_slope_days_per_year <- unname(coef(fit)["year"])
            out$year_p_value <- cf["year", 4]
          }
          if ("net_hours_100" %in% rownames(cf)) {
            out$net_hours_slope_days_per_100h <- unname(coef(fit)["net_hours_100"])
            out$net_hours_p_value <- cf["net_hours_100", 4]
          }
          if (nrow(cf) >= 4) {
            out$model_r_squared <- 1 - fit$deviance / fit$null.deviance
            out$model_aic <- AIC(fit)
          }
        }
      }
      glm_rows[[k]] <- out
    }
  }
}

glm_assoc <- bind_rows(glm_rows) %>%
  mutate(
    temp_effect_direction = ifelse(
      is.na(temp_slope_days_per_degC),
      NA_character_,
      ifelse(temp_slope_days_per_degC < 0, "warmer_month -> earlier_timing", "warmer_month -> later_timing")
    ),
    year_effect_direction = ifelse(
      is.na(year_slope_days_per_year),
      NA_character_,
      ifelse(year_slope_days_per_year < 0, "earlier_over_years", "later_over_years")
    ),
    net_hours_effect_direction = ifelse(
      is.na(net_hours_slope_days_per_100h),
      NA_character_,
      ifelse(net_hours_slope_days_per_100h < 0, "higher_effort -> earlier_timing", "higher_effort -> later_timing")
    )
  )

glm_assoc <- glm_assoc %>%
  mutate(
    temp_p_value_fdr_global = p.adjust(temp_p_value, method = "BH"),
    year_p_value_fdr_global = p.adjust(year_p_value, method = "BH"),
    net_hours_p_value_fdr_global = p.adjust(net_hours_p_value, method = "BH")
  ) %>%
  group_by(metric_label) %>%
  mutate(
    temp_p_value_fdr_within_metric = p.adjust(temp_p_value, method = "BH"),
    year_p_value_fdr_within_metric = p.adjust(year_p_value, method = "BH"),
    net_hours_p_value_fdr_within_metric = p.adjust(net_hours_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(temp_p_value)

write.csv(glm_assoc, out_glm_assoc, row.names = FALSE, na = "")

glm_summary <- glm_assoc %>%
  summarize(
    total_tests = sum(!is.na(temp_p_value)),
    temp_p_lt_0_05 = sum(temp_p_value < 0.05, na.rm = TRUE),
    temp_fdr_global_lt_0_10 = sum(temp_p_value_fdr_global < 0.10, na.rm = TRUE),
    temp_fdr_global_lt_0_05 = sum(temp_p_value_fdr_global < 0.05, na.rm = TRUE),
    year_p_lt_0_05 = sum(year_p_value < 0.05, na.rm = TRUE),
    year_fdr_global_lt_0_10 = sum(year_p_value_fdr_global < 0.10, na.rm = TRUE),
    year_fdr_global_lt_0_05 = sum(year_p_value_fdr_global < 0.05, na.rm = TRUE),
    net_hours_p_lt_0_05 = sum(net_hours_p_value < 0.05, na.rm = TRUE),
    net_hours_fdr_global_lt_0_10 = sum(net_hours_p_value_fdr_global < 0.10, na.rm = TRUE),
    net_hours_fdr_global_lt_0_05 = sum(net_hours_p_value_fdr_global < 0.05, na.rm = TRUE)
  )
write.csv(glm_summary, out_glm_summary, row.names = FALSE, na = "")

glm_sig <- glm_assoc %>%
  filter(!is.na(temp_p_value)) %>%
  filter(temp_p_value_fdr_global < 0.10 | temp_p_value < 0.05) %>%
  arrange(temp_p_value_fdr_global, temp_p_value)
write.csv(glm_sig, out_glm_sig, row.names = FALSE, na = "")

# 4) GAM sensitivity check, effort-adjusted: timing ~ year + net_hours_100 + s(temp)
# Run GAM for rows where GLM had valid fits.
glm_valid <- glm_assoc %>% filter(!is.na(temp_p_value))
gam_rows <- vector("list", nrow(glm_valid))

for (i in seq_len(nrow(glm_valid))) {
  r <- glm_valid[i, ]

  d_sp <- phen_long %>%
    filter(SPECIES == r$SPECIES, metric_label == r$metric_label) %>%
    select(year, timing_julian_day, net_hours_100)
  d_rg <- weather_wide %>% filter(region == r$region)
  merged <- inner_join(d_sp, d_rg, by = "year")

  mv <- names(month_labels)[match(r$temp_month, month_labels)]
  dd <- merged %>%
    transmute(
      year = year,
      timing_julian_day = timing_julian_day,
      net_hours_100 = net_hours_100,
      temp = .data[[mv]]
    ) %>%
    filter(!is.na(timing_julian_day), !is.na(temp), !is.na(net_hours_100), !is.na(year))

  out <- data.frame(
    SPECIES = r$SPECIES,
    SPNAME = r$SPNAME,
    phase_group = r$phase_group,
    metric_label = r$metric_label,
    region = r$region,
    temp_month = r$temp_month,
    n_years = nrow(dd),
    gam_temp_edf = NA_real_,
    gam_temp_p_value = NA_real_,
    gam_year_slope = NA_real_,
    gam_year_p_value = NA_real_,
    gam_net_hours_slope = NA_real_,
    gam_net_hours_p_value = NA_real_,
    gam_r_squared = NA_real_,
    gam_aic = NA_real_,
    glm_temp_slope = r$temp_slope_days_per_degC,
    glm_temp_p_value = r$temp_p_value,
    glm_aic = r$model_aic,
    delta_aic_gam_minus_glm = NA_real_,
    stringsAsFactors = FALSE
  )

  unique_temp <- length(unique(dd$temp))
  if (nrow(dd) >= 10 && stats::sd(dd$timing_julian_day) > 0 && stats::sd(dd$net_hours_100) > 0 && unique_temp >= 4) {
    k_temp <- max(3, min(5, unique_temp - 1))
    fit_gam <- safe_gam(timing_julian_day ~ year + net_hours_100 + s(temp, k = k_temp), dd)
    if (!is.null(fit_gam)) {
      sm <- summary(fit_gam)
      out$gam_aic <- AIC(fit_gam)
      out$delta_aic_gam_minus_glm <- out$gam_aic - out$glm_aic
      out$gam_r_squared <- sm$r.sq

      if (!is.null(sm$s.table) && nrow(sm$s.table) >= 1) {
        out$gam_temp_edf <- sm$s.table[1, "edf"]
        out$gam_temp_p_value <- sm$s.table[1, "p-value"]
      }
      if (!is.null(sm$p.table) && "year" %in% rownames(sm$p.table)) {
        out$gam_year_slope <- sm$p.table["year", "Estimate"]
        out$gam_year_p_value <- sm$p.table["year", "Pr(>|t|)"]
      }
      if (!is.null(sm$p.table) && "net_hours_100" %in% rownames(sm$p.table)) {
        out$gam_net_hours_slope <- sm$p.table["net_hours_100", "Estimate"]
        out$gam_net_hours_p_value <- sm$p.table["net_hours_100", "Pr(>|t|)"]
      }
    }
  }

  gam_rows[[i]] <- out
}

gam_assoc <- bind_rows(gam_rows) %>%
  mutate(
    gam_temp_p_value_fdr_global = p.adjust(gam_temp_p_value, method = "BH"),
    gam_year_p_value_fdr_global = p.adjust(gam_year_p_value, method = "BH"),
    gam_net_hours_p_value_fdr_global = p.adjust(gam_net_hours_p_value, method = "BH"),
    gam_temp_effect_direction = ifelse(
      is.na(glm_temp_slope),
      NA_character_,
      ifelse(glm_temp_slope < 0, "warmer_month -> earlier_timing", "warmer_month -> later_timing")
    )
  ) %>%
  arrange(gam_temp_p_value)

write.csv(gam_assoc, out_gam_assoc, row.names = FALSE, na = "")

gam_sig <- gam_assoc %>%
  filter(!is.na(gam_temp_p_value)) %>%
  filter(gam_temp_p_value_fdr_global < 0.10 | gam_temp_p_value < 0.05) %>%
  arrange(gam_temp_p_value_fdr_global, gam_temp_p_value)
write.csv(gam_sig, out_gam_sig, row.names = FALSE, na = "")

cat("Wrote:", out_weather_trend, "\n")
cat("Wrote:", out_phenology_trend, "\n")
cat("Wrote:", out_glm_assoc, "\n")
cat("Wrote:", out_glm_summary, "\n")
cat("Wrote:", out_glm_sig, "\n")
cat("Wrote:", out_gam_assoc, "\n")
cat("Wrote:", out_gam_sig, "\n")

cat("\nGLM association tests:", sum(!is.na(glm_assoc$temp_p_value)), "\n")
cat("GLM temp p<0.05:", sum(glm_assoc$temp_p_value < 0.05, na.rm = TRUE), "\n")
cat("GLM temp FDR<0.10:", sum(glm_assoc$temp_p_value_fdr_global < 0.10, na.rm = TRUE), "\n")
cat("GLM year p<0.05:", sum(glm_assoc$year_p_value < 0.05, na.rm = TRUE), "\n")
cat("GLM year FDR<0.10:", sum(glm_assoc$year_p_value_fdr_global < 0.10, na.rm = TRUE), "\n")
cat("GLM net_hours p<0.05:", sum(glm_assoc$net_hours_p_value < 0.05, na.rm = TRUE), "\n")
cat("GLM net_hours FDR<0.10:", sum(glm_assoc$net_hours_p_value_fdr_global < 0.10, na.rm = TRUE), "\n")

cat("\nGAM tests:", sum(!is.na(gam_assoc$gam_temp_p_value)), "\n")
cat("GAM temp smooth p<0.05:", sum(gam_assoc$gam_temp_p_value < 0.05, na.rm = TRUE), "\n")
cat("GAM temp smooth FDR<0.10:", sum(gam_assoc$gam_temp_p_value_fdr_global < 0.10, na.rm = TRUE), "\n")
cat("GAM net_hours p<0.05:", sum(gam_assoc$gam_net_hours_p_value < 0.05, na.rm = TRUE), "\n")
cat("GAM net_hours FDR<0.10:", sum(gam_assoc$gam_net_hours_p_value_fdr_global < 0.10, na.rm = TRUE), "\n")
