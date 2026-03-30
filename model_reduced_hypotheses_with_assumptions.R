options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(mgcv)
})

# Inputs
phenology_file <- "rpbo_focal_species_phenology_jul_oct.csv"
weather_wide_file <- "bc_weather_region_monthly_wide_jul_oct.csv"
net_hours_file <- "nethours.xlsx"

# Outputs
out_full_diag <- "model_full_glm_assumption_diagnostics.csv"
out_full_diag_summary <- "model_full_glm_assumption_summary.csv"
out_reduced <- "model_reduced_hypothesis_results.csv"
out_reduced_sig <- "model_reduced_hypothesis_significant.csv"
out_reduced_assump_summary <- "model_reduced_hypothesis_assumption_summary.csv"

set.seed(462)

metrics_lookup <- data.frame(
  metric_code = c("first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"),
  metric_label = c("first_detection", "p10", "p50", "p90", "last_detection"),
  stringsAsFactors = FALSE
)

month_vars <- c("july_mean_min_temp", "august_mean_min_temp", "september_mean_min_temp", "october_mean_min_temp")
month_labels <- c("july", "august", "september", "october")
names(month_labels) <- month_vars

is_valid_metric_for_phase <- function(phase_group, metric_label) {
  (phase_group == "departure" & metric_label != "first_detection") |
    (phase_group == "fall_arrival" & metric_label != "last_detection")
}

# Reduced hypothesis mapping (pre-specified):
# Departure species: p10/p50/p90/last only (no first detection; banding starts after presence begins).
# Fall-arrival species: first/p10/p50/p90 only (no last detection; banding ends before departure finishes).
reduced_map <- data.frame(
  phase_group = c(
    "departure", "departure", "departure", "departure",
    "fall_arrival", "fall_arrival", "fall_arrival", "fall_arrival"
  ),
  metric_label = c(
    "p10", "p50", "p90", "last_detection",
    "first_detection", "p10", "p50", "p90"
  ),
  temp_month = c(
    "july", "august", "september", "october",
    "september", "september", "september", "october"
  ),
  stringsAsFactors = FALSE
)

month_to_var <- setNames(names(month_labels), month_labels)

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

diag_from_lm <- function(fit, n) {
  res <- residuals(fit)
  fv <- fitted(fit)

  shapiro_p <- NA_real_
  if (length(res) >= 3 && length(res) <= 5000 && stats::sd(res) > 0) {
    shapiro_p <- tryCatch(stats::shapiro.test(res)$p.value, error = function(e) NA_real_)
  }

  homo_p <- NA_real_
  if (length(unique(fv)) > 1 && length(unique(abs(res))) > 1) {
    homo_p <- tryCatch(stats::cor.test(abs(res), fv, method = "spearman")$p.value, error = function(e) NA_real_)
  }

  lag1_resid_cor <- NA_real_
  if (length(res) >= 3 && stats::sd(res[-length(res)]) > 0 && stats::sd(res[-1]) > 0) {
    lag1_resid_cor <- suppressWarnings(stats::cor(res[-length(res)], res[-1]))
  }

  max_cooksd <- tryCatch(max(stats::cooks.distance(fit), na.rm = TRUE), error = function(e) NA_real_)
  cooksd_cutoff <- 4 / n

  pass_normal <- !is.na(shapiro_p) && shapiro_p > 0.05
  pass_homosced <- !is.na(homo_p) && homo_p > 0.05
  pass_autocorr <- !is.na(lag1_resid_cor) && abs(lag1_resid_cor) < 0.30
  pass_influence <- !is.na(max_cooksd) && max_cooksd < cooksd_cutoff
  pass_all <- pass_normal && pass_homosced && pass_autocorr && pass_influence

  data.frame(
    shapiro_p = shapiro_p,
    homoscedasticity_spearman_p = homo_p,
    lag1_residual_correlation = lag1_resid_cor,
    max_cooks_distance = max_cooksd,
    cooksd_cutoff = cooksd_cutoff,
    pass_normality = pass_normal,
    pass_homoscedasticity = pass_homosced,
    pass_autocorrelation = pass_autocorr,
    pass_influence = pass_influence,
    pass_all_assumptions = pass_all,
    stringsAsFactors = FALSE
  )
}

permute_temp_p <- function(dd, n_perm = 2000) {
  fit_obs <- safe_lm(timing_julian_day ~ year + net_hours_100 + temp, dd)
  if (is.null(fit_obs)) return(NA_real_)

  sm_obs <- summary(fit_obs)
  if (!("temp" %in% rownames(sm_obs$coefficients))) return(NA_real_)
  t_obs <- sm_obs$coefficients["temp", "t value"]
  if (is.na(t_obs)) return(NA_real_)

  t_perm <- rep(NA_real_, n_perm)
  for (b in seq_len(n_perm)) {
    dd$perm_temp <- sample(dd$temp, replace = FALSE)
    fit_b <- safe_lm(timing_julian_day ~ year + net_hours_100 + perm_temp, dd)
    if (is.null(fit_b)) next
    sm_b <- summary(fit_b)
    if ("perm_temp" %in% rownames(sm_b$coefficients)) {
      t_perm[b] <- sm_b$coefficients["perm_temp", "t value"]
    }
  }

  t_perm <- t_perm[!is.na(t_perm)]
  if (!length(t_perm)) return(NA_real_)
  mean(abs(t_perm) >= abs(t_obs))
}

# Load data
phen <- read.csv(phenology_file, stringsAsFactors = FALSE)
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

required_weather <- c("region", "year", month_vars)
missing_weather <- setdiff(required_weather, names(weather_wide))
if (length(missing_weather) > 0) {
  stop(paste("Missing weather columns:", paste(missing_weather, collapse = ", ")))
}

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
  filter(is_valid_metric_for_phase(phase_group, metric_label)) %>%
  left_join(net_hours, by = "year")

# Helper to build analysis frame
build_dd <- function(species, metric_label, region, temp_var) {
  d_sp <- phen_long %>%
    filter(SPECIES == species, metric_label == !!metric_label) %>%
    select(year, timing_julian_day, net_hours_100)
  d_rg <- weather_wide %>% filter(region == !!region)
  merged <- inner_join(d_sp, d_rg, by = "year")
  merged %>%
    transmute(
      year = year,
      timing_julian_day = timing_julian_day,
      net_hours_100 = net_hours_100,
      temp = .data[[temp_var]]
    ) %>%
    filter(!is.na(year), !is.na(timing_julian_day), !is.na(net_hours_100), !is.na(temp))
}

# ------------------------------------------------------------------------------
# 1) Assumption diagnostics across the full model set (all 4 months x metrics).
# ------------------------------------------------------------------------------
keys_full <- phen_long %>%
  distinct(SPECIES, SPNAME, phase_group) %>%
  tidyr::crossing(
    metric_label = metrics_lookup$metric_label,
    region = sort(unique(weather_wide$region)),
    temp_var = month_vars
  ) %>%
  filter(is_valid_metric_for_phase(phase_group, metric_label))

full_rows <- vector("list", nrow(keys_full))
for (i in seq_len(nrow(keys_full))) {
  sp <- keys_full$SPECIES[i]
  mt <- keys_full$metric_label[i]
  rg <- keys_full$region[i]
  tv <- keys_full$temp_var[i]

  spname <- keys_full$SPNAME[i]
  phase <- keys_full$phase_group[i]
  dd <- build_dd(sp, mt, rg, tv)

  out <- data.frame(
    SPECIES = sp,
    SPNAME = spname,
    phase_group = phase,
    metric_label = mt,
    region = rg,
    temp_month = month_labels[[tv]],
    n_years = nrow(dd),
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 10 && stats::sd(dd$timing_julian_day) > 0 && stats::sd(dd$net_hours_100) > 0 && stats::sd(dd$temp) > 0) {
    fit <- safe_lm(timing_julian_day ~ year + net_hours_100 + temp, dd)
    if (!is.null(fit)) {
      out <- cbind(out, diag_from_lm(fit, nrow(dd)))
    } else {
      out <- cbind(out, data.frame(
        shapiro_p = NA_real_,
        homoscedasticity_spearman_p = NA_real_,
        lag1_residual_correlation = NA_real_,
        max_cooks_distance = NA_real_,
        cooksd_cutoff = NA_real_,
        pass_normality = FALSE,
        pass_homoscedasticity = FALSE,
        pass_autocorrelation = FALSE,
        pass_influence = FALSE,
        pass_all_assumptions = FALSE
      ))
    }
  } else {
    out <- cbind(out, data.frame(
      shapiro_p = NA_real_,
      homoscedasticity_spearman_p = NA_real_,
      lag1_residual_correlation = NA_real_,
      max_cooks_distance = NA_real_,
      cooksd_cutoff = NA_real_,
      pass_normality = FALSE,
      pass_homoscedasticity = FALSE,
      pass_autocorrelation = FALSE,
      pass_influence = FALSE,
      pass_all_assumptions = FALSE
    ))
  }
  full_rows[[i]] <- out
}

full_diag <- bind_rows(full_rows)
write.csv(full_diag, out_full_diag, row.names = FALSE, na = "")

full_summary <- data.frame(
  total_models = nrow(full_diag),
  models_with_enough_data = sum(full_diag$n_years >= 10),
  pass_normality = sum(full_diag$pass_normality, na.rm = TRUE),
  pass_homoscedasticity = sum(full_diag$pass_homoscedasticity, na.rm = TRUE),
  pass_autocorrelation = sum(full_diag$pass_autocorrelation, na.rm = TRUE),
  pass_influence = sum(full_diag$pass_influence, na.rm = TRUE),
  pass_all_assumptions = sum(full_diag$pass_all_assumptions, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(full_summary, out_full_diag_summary, row.names = FALSE, na = "")

# ------------------------------------------------------------------------------
# 2) Reduced hypothesis model set with diagnostics and robust p-values.
# ------------------------------------------------------------------------------
keys_reduced <- phen_long %>%
  distinct(SPECIES, SPNAME, phase_group, metric_label) %>%
  inner_join(reduced_map, by = c("phase_group", "metric_label")) %>%
  tidyr::crossing(region = sort(unique(weather_wide$region))) %>%
  mutate(temp_var = month_to_var[temp_month])

red_rows <- vector("list", nrow(keys_reduced))
for (i in seq_len(nrow(keys_reduced))) {
  r <- keys_reduced[i, ]
  dd <- build_dd(r$SPECIES, r$metric_label, r$region, r$temp_var)

  out <- data.frame(
    SPECIES = r$SPECIES,
    SPNAME = r$SPNAME,
    phase_group = r$phase_group,
    metric_label = r$metric_label,
    region = r$region,
    temp_month = r$temp_month,
    n_years = nrow(dd),
    temp_slope_days_per_degC = NA_real_,
    temp_p_value_glm = NA_real_,
    temp_p_value_permutation = NA_real_,
    year_slope_days_per_year = NA_real_,
    year_p_value_glm = NA_real_,
    net_hours_slope_days_per_100h = NA_real_,
    net_hours_p_value_glm = NA_real_,
    model_r_squared = NA_real_,
    model_aic = NA_real_,
    gam_temp_edf = NA_real_,
    gam_temp_p_value = NA_real_,
    gam_year_slope = NA_real_,
    gam_year_p_value = NA_real_,
    gam_net_hours_slope = NA_real_,
    gam_net_hours_p_value = NA_real_,
    gam_r_squared = NA_real_,
    gam_aic = NA_real_,
    delta_aic_gam_minus_glm = NA_real_,
    stringsAsFactors = FALSE
  )

  diag_blank <- data.frame(
    shapiro_p = NA_real_,
    homoscedasticity_spearman_p = NA_real_,
    lag1_residual_correlation = NA_real_,
    max_cooks_distance = NA_real_,
    cooksd_cutoff = NA_real_,
    pass_normality = FALSE,
    pass_homoscedasticity = FALSE,
    pass_autocorrelation = FALSE,
    pass_influence = FALSE,
    pass_all_assumptions = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(dd) >= 10 && stats::sd(dd$timing_julian_day) > 0 && stats::sd(dd$net_hours_100) > 0 && stats::sd(dd$temp) > 0) {
    fit_glm <- safe_glm(timing_julian_day ~ year + net_hours_100 + temp, dd)
    if (!is.null(fit_glm)) {
      sm <- summary(fit_glm)
      cf <- sm$coefficients
      if ("temp" %in% rownames(cf)) {
        out$temp_slope_days_per_degC <- unname(coef(fit_glm)["temp"])
        out$temp_p_value_glm <- cf["temp", 4]
      }
      if ("year" %in% rownames(cf)) {
        out$year_slope_days_per_year <- unname(coef(fit_glm)["year"])
        out$year_p_value_glm <- cf["year", 4]
      }
      if ("net_hours_100" %in% rownames(cf)) {
        out$net_hours_slope_days_per_100h <- unname(coef(fit_glm)["net_hours_100"])
        out$net_hours_p_value_glm <- cf["net_hours_100", 4]
      }
      out$model_r_squared <- 1 - fit_glm$deviance / fit_glm$null.deviance
      out$model_aic <- AIC(fit_glm)

      fit_lm_for_diag <- safe_lm(timing_julian_day ~ year + net_hours_100 + temp, dd)
      if (!is.null(fit_lm_for_diag)) {
        diag_blank <- diag_from_lm(fit_lm_for_diag, nrow(dd))
      }

      out$temp_p_value_permutation <- permute_temp_p(dd, n_perm = 2000)
    }

    # GAM sensitivity check: allow non-linear temp effect.
    unique_temp <- length(unique(dd$temp))
    if (unique_temp >= 4) {
      k_temp <- max(3, min(5, unique_temp - 1))
      fit_gam <- safe_gam(timing_julian_day ~ year + net_hours_100 + s(temp, k = k_temp), dd)
      if (!is.null(fit_gam)) {
        smg <- summary(fit_gam)
        out$gam_aic <- AIC(fit_gam)
        out$delta_aic_gam_minus_glm <- out$gam_aic - out$model_aic
        out$gam_r_squared <- smg$r.sq
        if (!is.null(smg$s.table) && nrow(smg$s.table) >= 1) {
          out$gam_temp_edf <- smg$s.table[1, "edf"]
          out$gam_temp_p_value <- smg$s.table[1, "p-value"]
        }
        if (!is.null(smg$p.table) && "year" %in% rownames(smg$p.table)) {
          out$gam_year_slope <- smg$p.table["year", "Estimate"]
          out$gam_year_p_value <- smg$p.table["year", "Pr(>|t|)"]
        }
        if (!is.null(smg$p.table) && "net_hours_100" %in% rownames(smg$p.table)) {
          out$gam_net_hours_slope <- smg$p.table["net_hours_100", "Estimate"]
          out$gam_net_hours_p_value <- smg$p.table["net_hours_100", "Pr(>|t|)"]
        }
      }
    }
  }

  red_rows[[i]] <- cbind(out, diag_blank)
}

reduced <- bind_rows(red_rows) %>%
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
  ) %>%
  mutate(
    temp_p_fdr_global_glm = p.adjust(temp_p_value_glm, method = "BH"),
    temp_p_fdr_global_permutation = p.adjust(temp_p_value_permutation, method = "BH"),
    temp_p_fdr_global_gam = p.adjust(gam_temp_p_value, method = "BH"),
    year_p_fdr_global_glm = p.adjust(year_p_value_glm, method = "BH"),
    net_hours_p_fdr_global_glm = p.adjust(net_hours_p_value_glm, method = "BH"),
    gam_net_hours_p_fdr_global = p.adjust(gam_net_hours_p_value, method = "BH")
  ) %>%
  group_by(metric_label) %>%
  mutate(
    temp_p_fdr_within_metric_glm = p.adjust(temp_p_value_glm, method = "BH"),
    temp_p_fdr_within_metric_permutation = p.adjust(temp_p_value_permutation, method = "BH"),
    temp_p_fdr_within_metric_gam = p.adjust(gam_temp_p_value, method = "BH"),
    year_p_fdr_within_metric_glm = p.adjust(year_p_value_glm, method = "BH"),
    net_hours_p_fdr_within_metric_glm = p.adjust(net_hours_p_value_glm, method = "BH"),
    gam_net_hours_p_fdr_within_metric = p.adjust(gam_net_hours_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(temp_p_value_glm)

write.csv(reduced, out_reduced, row.names = FALSE, na = "")

reduced_sig <- reduced %>%
  filter(!is.na(temp_p_value_glm)) %>%
  filter(
    temp_p_fdr_global_glm < 0.10 |
      temp_p_fdr_global_permutation < 0.10 |
      temp_p_fdr_global_gam < 0.10 |
      temp_p_value_glm < 0.05 |
      temp_p_value_permutation < 0.05 |
      gam_temp_p_value < 0.05
  ) %>%
  arrange(temp_p_fdr_global_glm, temp_p_fdr_global_permutation, temp_p_fdr_global_gam, temp_p_value_glm)

write.csv(reduced_sig, out_reduced_sig, row.names = FALSE, na = "")

reduced_assump_summary <- data.frame(
  total_reduced_models = nrow(reduced),
  pass_normality = sum(reduced$pass_normality, na.rm = TRUE),
  pass_homoscedasticity = sum(reduced$pass_homoscedasticity, na.rm = TRUE),
  pass_autocorrelation = sum(reduced$pass_autocorrelation, na.rm = TRUE),
  pass_influence = sum(reduced$pass_influence, na.rm = TRUE),
  pass_all_assumptions = sum(reduced$pass_all_assumptions, na.rm = TRUE),
  glm_temp_p_lt_0_05 = sum(reduced$temp_p_value_glm < 0.05, na.rm = TRUE),
  glm_temp_fdr_global_lt_0_10 = sum(reduced$temp_p_fdr_global_glm < 0.10, na.rm = TRUE),
  glm_net_hours_p_lt_0_05 = sum(reduced$net_hours_p_value_glm < 0.05, na.rm = TRUE),
  glm_net_hours_fdr_global_lt_0_10 = sum(reduced$net_hours_p_fdr_global_glm < 0.10, na.rm = TRUE),
  perm_temp_p_lt_0_05 = sum(reduced$temp_p_value_permutation < 0.05, na.rm = TRUE),
  perm_temp_fdr_global_lt_0_10 = sum(reduced$temp_p_fdr_global_permutation < 0.10, na.rm = TRUE),
  gam_temp_p_lt_0_05 = sum(reduced$gam_temp_p_value < 0.05, na.rm = TRUE),
  gam_temp_fdr_global_lt_0_10 = sum(reduced$temp_p_fdr_global_gam < 0.10, na.rm = TRUE),
  gam_net_hours_p_lt_0_05 = sum(reduced$gam_net_hours_p_value < 0.05, na.rm = TRUE),
  gam_net_hours_fdr_global_lt_0_10 = sum(reduced$gam_net_hours_p_fdr_global < 0.10, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(reduced_assump_summary, out_reduced_assump_summary, row.names = FALSE, na = "")

cat("Wrote:", out_full_diag, "\n")
cat("Wrote:", out_full_diag_summary, "\n")
cat("Wrote:", out_reduced, "\n")
cat("Wrote:", out_reduced_sig, "\n")
cat("Wrote:", out_reduced_assump_summary, "\n")

cat("\nFull model set assumption pass (all):", full_summary$pass_all_assumptions, "/", full_summary$total_models, "\n")
cat("Reduced model set assumption pass (all):", reduced_assump_summary$pass_all_assumptions, "/", reduced_assump_summary$total_reduced_models, "\n")
cat("Reduced GLM temp p<0.05:", reduced_assump_summary$glm_temp_p_lt_0_05, "\n")
cat("Reduced GLM temp FDR<0.10:", reduced_assump_summary$glm_temp_fdr_global_lt_0_10, "\n")
cat("Reduced GLM net_hours p<0.05:", reduced_assump_summary$glm_net_hours_p_lt_0_05, "\n")
cat("Reduced GLM net_hours FDR<0.10:", reduced_assump_summary$glm_net_hours_fdr_global_lt_0_10, "\n")
cat("Reduced permutation temp p<0.05:", reduced_assump_summary$perm_temp_p_lt_0_05, "\n")
cat("Reduced permutation temp FDR<0.10:", reduced_assump_summary$perm_temp_fdr_global_lt_0_10, "\n")
cat("Reduced GAM temp p<0.05:", reduced_assump_summary$gam_temp_p_lt_0_05, "\n")
cat("Reduced GAM temp FDR<0.10:", reduced_assump_summary$gam_temp_fdr_global_lt_0_10, "\n")
cat("Reduced GAM net_hours p<0.05:", reduced_assump_summary$gam_net_hours_p_lt_0_05, "\n")
cat("Reduced GAM net_hours FDR<0.10:", reduced_assump_summary$gam_net_hours_fdr_global_lt_0_10, "\n")
