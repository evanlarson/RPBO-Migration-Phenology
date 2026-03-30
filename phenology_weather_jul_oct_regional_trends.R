options(stringsAsFactors = FALSE)

# Inputs
rpbo_input <- "rpbo_clean_standardized.csv"
weather_pattern <- "_temp_clean\\.csv$"

# Outputs
species_group_output <- "focal_species_groups_jul_oct.csv"
rpbo_filtered_output <- "rpbo_jul_oct_only.csv"
phenology_output <- "rpbo_focal_species_phenology_jul_oct.csv"
weather_station_monthly_output <- "bc_weather_station_monthly_jul_oct.csv"
weather_region_monthly_output <- "bc_weather_region_monthly_jul_oct.csv"
weather_region_wide_output <- "bc_weather_region_monthly_wide_jul_oct.csv"
phenology_year_trend_output <- "phenology_year_trends_jul_oct.csv"
weather_year_trend_output <- "regional_weather_year_trends_jul_oct.csv"
phenology_temp_monthly_models_output <- "phenology_vs_monthly_temp_models_jul_oct.csv"
phenology_temp_full_models_output <- "phenology_vs_all_months_models_jul_oct.csv"

# Region mapping requested by user.
region_map <- data.frame(
  station = c(
    "victoria", "vancouver", "courtenay",
    "porthardy", "bellacoola", "princerupert",
    "kamloops", "princeton", "pemberton",
    "princegeorge", "williamslake", "smithers"
  ),
  region = c(
    "southern_coastal_bc", "southern_coastal_bc", "southern_coastal_bc",
    "central_coast", "central_coast", "central_coast",
    "south_interior", "south_interior", "south_interior",
    "central_interior_sub_boreal", "central_interior_sub_boreal", "central_interior_sub_boreal"
  ),
  stringsAsFactors = FALSE
)

# Helpers
weighted_quantile <- function(x, w, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  out <- sapply(probs, function(p) x[which(cw >= p)[1]])
  names(out) <- paste0("q", probs * 100)
  out
}

safe_fit <- function(formula, data) {
  fit <- try(lm(formula, data = data), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

month_name <- function(m) {
  ifelse(m == 7, "july",
    ifelse(m == 8, "august",
      ifelse(m == 9, "september",
        ifelse(m == 10, "october", NA_character_)
      )
    )
  )
}

# 1) Load RPBO and keep Jul-Oct only.
rpbo <- read.csv(rpbo_input, check.names = FALSE, stringsAsFactors = FALSE)
need_rpbo <- c("SPECIES", "SPNAME", "Banded", "year_std", "julian_day", "date_iso")
miss_rpbo <- setdiff(need_rpbo, names(rpbo))
if (length(miss_rpbo) > 0) {
  stop(paste("Missing RPBO columns:", paste(miss_rpbo, collapse = ", ")))
}

rpbo$date <- as.Date(rpbo$date_iso)
rpbo$month <- as.integer(format(rpbo$date, "%m"))

rpbo_jul_oct <- rpbo[
  !is.na(rpbo$month) &
    rpbo$month %in% c(7, 8, 9, 10) &
    !is.na(rpbo$Banded) &
    rpbo$Banded > 0 &
    !is.na(rpbo$year_std) &
    !is.na(rpbo$julian_day),
]
if (nrow(rpbo_jul_oct) == 0) {
  stop("No Jul-Oct RPBO rows with Banded > 0 were found.")
}

write.csv(rpbo_jul_oct, rpbo_filtered_output, row.names = FALSE, na = "")

# 2) Define 3 departure + 3 fall-arrival species from Jul-Oct timing profile.
species_split <- split(rpbo_jul_oct, rpbo_jul_oct$SPECIES)
species_rows <- lapply(species_split, function(d) {
  q <- weighted_quantile(d$julian_day, d$Banded, probs = c(0.1, 0.5, 0.9))
  data.frame(
    SPECIES = d$SPECIES[1],
    SPNAME = d$SPNAME[1],
    total_banded = sum(d$Banded),
    years_with_data = length(unique(d$year_std)),
    overall_q10_julian = as.integer(q["q10"]),
    overall_q50_julian = as.integer(q["q50"]),
    overall_q90_julian = as.integer(q["q90"]),
    stringsAsFactors = FALSE
  )
})
species_summary <- do.call(rbind, species_rows)
species_summary <- species_summary[order(species_summary$overall_q50_julian), ]

n_species <- nrow(species_summary)
if (n_species != 6) {
  warning(sprintf("Expected 6 focal species but found %d in Jul-Oct filtered data.", n_species))
}
split_n <- floor(n_species / 2)
departure_species <- species_summary$SPECIES[seq_len(split_n)]
arrival_species <- species_summary$SPECIES[(split_n + 1):n_species]

species_summary$phase_group <- ifelse(
  species_summary$SPECIES %in% arrival_species,
  "fall_arrival",
  "departure"
)
write.csv(species_summary, species_group_output, row.names = FALSE, na = "")

# 3) Per-species per-year phenology metrics (Jul-Oct only).
key <- paste(rpbo_jul_oct$SPECIES, rpbo_jul_oct$year_std, sep = "__")
grouped <- split(rpbo_jul_oct, key)

phenology_rows <- lapply(grouped, function(d) {
  q <- weighted_quantile(d$julian_day, d$Banded, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  species <- d$SPECIES[1]
  phase <- if (species %in% arrival_species) "fall_arrival" else "departure"
  phase_metric <- if (phase == "fall_arrival") "q10_julian" else "q90_julian"
  phase_julian <- if (phase == "fall_arrival") as.integer(q["q10"]) else as.integer(q["q90"])

  data.frame(
    year = as.integer(d$year_std[1]),
    SPECIES = species,
    SPNAME = d$SPNAME[1],
    phase_group = phase,
    phase_metric = phase_metric,
    phase_julian_day = phase_julian,
    first_julian_day = min(d$julian_day),
    q10_julian_day = as.integer(q["q10"]),
    q25_julian_day = as.integer(q["q25"]),
    q50_julian_day = as.integer(q["q50"]),
    q75_julian_day = as.integer(q["q75"]),
    q90_julian_day = as.integer(q["q90"]),
    last_julian_day = max(d$julian_day),
    total_banded = sum(d$Banded),
    sampled_days = nrow(d),
    stringsAsFactors = FALSE
  )
})
phenology <- do.call(rbind, phenology_rows)
phenology <- phenology[order(phenology$SPECIES, phenology$year), ]
write.csv(phenology, phenology_output, row.names = FALSE, na = "")

analysis_years <- sort(unique(phenology$year))

# 4) Weather monthly means by station (Jul-Oct only), then by requested region.
weather_files <- list.files(pattern = weather_pattern)
if (length(weather_files) == 0) {
  stop("No *_temp_clean.csv files found.")
}

station_monthly_rows <- list()
for (f in weather_files) {
  station <- sub("_temp_clean\\.csv$", "", f)
  reg_idx <- match(station, region_map$station)
  if (is.na(reg_idx)) {
    next
  }
  region <- region_map$region[reg_idx]

  w <- read.csv(f, stringsAsFactors = FALSE)
  need_w <- c("date", "year", "min_temperature")
  miss_w <- setdiff(need_w, names(w))
  if (length(miss_w) > 0) {
    stop(paste("Missing weather columns in", f, ":", paste(miss_w, collapse = ", ")))
  }

  w$date <- as.Date(w$date)
  w$month <- as.integer(format(w$date, "%m"))
  w <- w[
    !is.na(w$year) &
      w$year %in% analysis_years &
      !is.na(w$month) &
      w$month %in% c(7, 8, 9, 10) &
      !is.na(w$min_temperature),
  ]

  if (!nrow(w)) next

  keys <- paste(w$year, w$month, sep = "__")
  by_ym <- split(w, keys)

  rows <- lapply(by_ym, function(d) {
    data.frame(
      station = station,
      region = region,
      year = as.integer(d$year[1]),
      month = as.integer(d$month[1]),
      month_name = month_name(as.integer(d$month[1])),
      n_days = nrow(d),
      mean_min_temp = mean(d$min_temperature, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  station_monthly_rows <- c(station_monthly_rows, rows)
}

weather_station_monthly <- do.call(rbind, station_monthly_rows)
weather_station_monthly <- weather_station_monthly[order(
  weather_station_monthly$region, weather_station_monthly$station,
  weather_station_monthly$year, weather_station_monthly$month
), ]
write.csv(weather_station_monthly, weather_station_monthly_output, row.names = FALSE, na = "")

rm_keys <- paste(weather_station_monthly$region, weather_station_monthly$year, weather_station_monthly$month, sep = "__")
region_groups <- split(weather_station_monthly, rm_keys)
region_monthly_rows <- lapply(region_groups, function(d) {
  data.frame(
    region = d$region[1],
    year = as.integer(d$year[1]),
    month = as.integer(d$month[1]),
    month_name = d$month_name[1],
    n_station_values = nrow(d),
    region_mean_min_temp = mean(d$mean_min_temp, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

weather_region_monthly <- do.call(rbind, region_monthly_rows)
weather_region_monthly <- weather_region_monthly[order(
  weather_region_monthly$region, weather_region_monthly$year, weather_region_monthly$month
), ]
write.csv(weather_region_monthly, weather_region_monthly_output, row.names = FALSE, na = "")

# Region-year wide table with Jul/Aug/Sep/Oct bins.
ry_keys <- paste(weather_region_monthly$region, weather_region_monthly$year, sep = "__")
ry_groups <- split(weather_region_monthly, ry_keys)
weather_region_wide_rows <- lapply(ry_groups, function(d) {
  out <- data.frame(
    region = d$region[1],
    year = as.integer(d$year[1]),
    july_mean_min_temp = NA_real_,
    august_mean_min_temp = NA_real_,
    september_mean_min_temp = NA_real_,
    october_mean_min_temp = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(d))) {
    m <- d$month[i]
    if (m == 7) out$july_mean_min_temp <- d$region_mean_min_temp[i]
    if (m == 8) out$august_mean_min_temp <- d$region_mean_min_temp[i]
    if (m == 9) out$september_mean_min_temp <- d$region_mean_min_temp[i]
    if (m == 10) out$october_mean_min_temp <- d$region_mean_min_temp[i]
  }
  out
})
weather_region_wide <- do.call(rbind, weather_region_wide_rows)
weather_region_wide <- weather_region_wide[order(weather_region_wide$region, weather_region_wide$year), ]
write.csv(weather_region_wide, weather_region_wide_output, row.names = FALSE, na = "")

# 5) Are birds shifting timing over years? (phase_julian_day ~ year)
sp_groups <- split(phenology, phenology$SPECIES)
phenology_trend_rows <- lapply(sp_groups, function(d) {
  keep <- !is.na(d$phase_julian_day) & !is.na(d$year)
  d <- d[keep, ]
  fit <- safe_fit(phase_julian_day ~ year, d)
  if (is.null(fit) || nrow(d) < 8 || stats::sd(d$phase_julian_day) == 0) {
    return(data.frame(
      SPECIES = d$SPECIES[1],
      SPNAME = d$SPNAME[1],
      phase_group = d$phase_group[1],
      n_years = nrow(d),
      slope_days_per_year = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  sm <- summary(fit)
  data.frame(
    SPECIES = d$SPECIES[1],
    SPNAME = d$SPNAME[1],
    phase_group = d$phase_group[1],
    n_years = nrow(d),
    slope_days_per_year = unname(coef(fit)[2]),
    p_value = sm$coefficients[2, 4],
    r_squared = sm$r.squared,
    stringsAsFactors = FALSE
  )
})
phenology_year_trends <- do.call(rbind, phenology_trend_rows)
phenology_year_trends <- phenology_year_trends[order(phenology_year_trends$p_value), ]
phenology_year_trends$p_value_fdr <- p.adjust(phenology_year_trends$p_value, method = "BH")
write.csv(phenology_year_trends, phenology_year_trend_output, row.names = FALSE, na = "")

# 6) Is climate warming in each region/month? (temp ~ year)
rw_groups <- split(weather_region_monthly, paste(weather_region_monthly$region, weather_region_monthly$month_name, sep = "__"))
weather_trend_rows <- lapply(rw_groups, function(d) {
  keep <- !is.na(d$region_mean_min_temp) & !is.na(d$year)
  d <- d[keep, ]
  fit <- safe_fit(region_mean_min_temp ~ year, d)
  if (is.null(fit) || nrow(d) < 8 || stats::sd(d$region_mean_min_temp) == 0) {
    return(data.frame(
      region = d$region[1],
      month_name = d$month_name[1],
      n_years = nrow(d),
      slope_degC_per_year = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  sm <- summary(fit)
  data.frame(
    region = d$region[1],
    month_name = d$month_name[1],
    n_years = nrow(d),
    slope_degC_per_year = unname(coef(fit)[2]),
    p_value = sm$coefficients[2, 4],
    r_squared = sm$r.squared,
    stringsAsFactors = FALSE
  )
})
weather_year_trends <- do.call(rbind, weather_trend_rows)
weather_year_trends <- weather_year_trends[order(weather_year_trends$p_value), ]
weather_year_trends$p_value_fdr <- p.adjust(weather_year_trends$p_value, method = "BH")
write.csv(weather_year_trends, weather_year_trend_output, row.names = FALSE, na = "")

# 7) Month-by-month temperature effects on timing, controlling for year.
months_interest <- c("july", "august", "september", "october")
regions <- sort(unique(weather_region_wide$region))
species_list <- sort(unique(phenology$SPECIES))

month_model_rows <- list()
full_model_rows <- list()

for (sp in species_list) {
  d_sp <- phenology[phenology$SPECIES == sp, ]
  for (rg in regions) {
    d_rg <- weather_region_wide[weather_region_wide$region == rg, ]
    merged <- merge(d_sp, d_rg, by = "year", all = FALSE)
    merged <- merged[order(merged$year), ]

    # Month-by-month models: phase ~ year + month_temp
    for (m in months_interest) {
      temp_col <- paste0(m, "_mean_min_temp")
      keep <- !is.na(merged$phase_julian_day) & !is.na(merged$year) & !is.na(merged[[temp_col]])
      dd <- merged[keep, ]
      if (nrow(dd) < 10 || stats::sd(dd$phase_julian_day) == 0 || stats::sd(dd[[temp_col]]) == 0) {
        month_model_rows[[length(month_model_rows) + 1]] <- data.frame(
          SPECIES = sp,
          SPNAME = d_sp$SPNAME[1],
          phase_group = d_sp$phase_group[1],
          region = rg,
          month_name = m,
          n_years = nrow(dd),
          temp_slope_days_per_degC = NA_real_,
          temp_p_value = NA_real_,
          year_slope_days_per_year = NA_real_,
          year_p_value = NA_real_,
          model_r_squared = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      fit <- lm(phase_julian_day ~ year + dd[[temp_col]], data = dd)
      sm <- summary(fit)
      coefs <- sm$coefficients
      month_model_rows[[length(month_model_rows) + 1]] <- data.frame(
        SPECIES = sp,
        SPNAME = d_sp$SPNAME[1],
        phase_group = d_sp$phase_group[1],
        region = rg,
        month_name = m,
        n_years = nrow(dd),
        temp_slope_days_per_degC = unname(coef(fit)[3]),
        temp_p_value = coefs[3, 4],
        year_slope_days_per_year = unname(coef(fit)[2]),
        year_p_value = coefs[2, 4],
        model_r_squared = sm$r.squared,
        stringsAsFactors = FALSE
      )
    }

    # Full 4-month model: phase ~ year + Jul + Aug + Sep + Oct
    needed <- c(
      "phase_julian_day", "year",
      "july_mean_min_temp", "august_mean_min_temp", "september_mean_min_temp", "october_mean_min_temp"
    )
    keep_full <- complete.cases(merged[, needed])
    dfull <- merged[keep_full, ]
    if (nrow(dfull) < 12 || stats::sd(dfull$phase_julian_day) == 0) {
      full_model_rows[[length(full_model_rows) + 1]] <- data.frame(
        SPECIES = sp,
        SPNAME = d_sp$SPNAME[1],
        phase_group = d_sp$phase_group[1],
        region = rg,
        n_years = nrow(dfull),
        year_slope_days_per_year = NA_real_,
        year_p_value = NA_real_,
        july_slope_days_per_degC = NA_real_,
        july_p_value = NA_real_,
        august_slope_days_per_degC = NA_real_,
        august_p_value = NA_real_,
        september_slope_days_per_degC = NA_real_,
        september_p_value = NA_real_,
        october_slope_days_per_degC = NA_real_,
        october_p_value = NA_real_,
        model_r_squared = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    fit_full <- safe_fit(
      phase_julian_day ~ year + july_mean_min_temp + august_mean_min_temp + september_mean_min_temp + october_mean_min_temp,
      dfull
    )
    if (is.null(fit_full)) {
      full_model_rows[[length(full_model_rows) + 1]] <- data.frame(
        SPECIES = sp,
        SPNAME = d_sp$SPNAME[1],
        phase_group = d_sp$phase_group[1],
        region = rg,
        n_years = nrow(dfull),
        year_slope_days_per_year = NA_real_,
        year_p_value = NA_real_,
        july_slope_days_per_degC = NA_real_,
        july_p_value = NA_real_,
        august_slope_days_per_degC = NA_real_,
        august_p_value = NA_real_,
        september_slope_days_per_degC = NA_real_,
        september_p_value = NA_real_,
        october_slope_days_per_degC = NA_real_,
        october_p_value = NA_real_,
        model_r_squared = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      smf <- summary(fit_full)
      cf <- smf$coefficients
      full_model_rows[[length(full_model_rows) + 1]] <- data.frame(
        SPECIES = sp,
        SPNAME = d_sp$SPNAME[1],
        phase_group = d_sp$phase_group[1],
        region = rg,
        n_years = nrow(dfull),
        year_slope_days_per_year = unname(coef(fit_full)["year"]),
        year_p_value = cf["year", 4],
        july_slope_days_per_degC = unname(coef(fit_full)["july_mean_min_temp"]),
        july_p_value = cf["july_mean_min_temp", 4],
        august_slope_days_per_degC = unname(coef(fit_full)["august_mean_min_temp"]),
        august_p_value = cf["august_mean_min_temp", 4],
        september_slope_days_per_degC = unname(coef(fit_full)["september_mean_min_temp"]),
        september_p_value = cf["september_mean_min_temp", 4],
        october_slope_days_per_degC = unname(coef(fit_full)["october_mean_min_temp"]),
        october_p_value = cf["october_mean_min_temp", 4],
        model_r_squared = smf$r.squared,
        stringsAsFactors = FALSE
      )
    }
  }
}

month_models <- do.call(rbind, month_model_rows)
month_models <- month_models[order(month_models$temp_p_value), ]
month_models$temp_p_value_fdr <- p.adjust(month_models$temp_p_value, method = "BH")
write.csv(month_models, phenology_temp_monthly_models_output, row.names = FALSE, na = "")

full_models <- do.call(rbind, full_model_rows)
full_models <- full_models[order(full_models$year_p_value), ]

full_models$july_p_value_fdr <- p.adjust(full_models$july_p_value, method = "BH")
full_models$august_p_value_fdr <- p.adjust(full_models$august_p_value, method = "BH")
full_models$september_p_value_fdr <- p.adjust(full_models$september_p_value, method = "BH")
full_models$october_p_value_fdr <- p.adjust(full_models$october_p_value, method = "BH")
full_models$year_p_value_fdr <- p.adjust(full_models$year_p_value, method = "BH")
write.csv(full_models, phenology_temp_full_models_output, row.names = FALSE, na = "")

cat("Wrote:", species_group_output, "\n")
cat("Wrote:", rpbo_filtered_output, "\n")
cat("Wrote:", phenology_output, "\n")
cat("Wrote:", weather_station_monthly_output, "\n")
cat("Wrote:", weather_region_monthly_output, "\n")
cat("Wrote:", weather_region_wide_output, "\n")
cat("Wrote:", phenology_year_trend_output, "\n")
cat("Wrote:", weather_year_trend_output, "\n")
cat("Wrote:", phenology_temp_monthly_models_output, "\n")
cat("Wrote:", phenology_temp_full_models_output, "\n")

cat("Years analyzed:", min(analysis_years), "-", max(analysis_years), "(n =", length(analysis_years), ")\n")
cat("Species analyzed:", n_species, "\n")
cat("Jul-Oct RPBO rows:", nrow(rpbo_jul_oct), "\n")
cat("Phenology rows:", nrow(phenology), "\n")
cat("Station-month weather rows:", nrow(weather_station_monthly), "\n")
cat("Region-month weather rows:", nrow(weather_region_monthly), "\n")
cat("Month-by-month models:", nrow(month_models), "\n")
cat("Full 4-month models:", nrow(full_models), "\n")
