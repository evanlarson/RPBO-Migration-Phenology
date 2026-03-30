options(stringsAsFactors = FALSE)

rpbo_input <- "rpbo_clean_standardized.csv"
species_group_output <- "focal_species_groups.csv"
phenology_output <- "rpbo_focal_species_phenology.csv"
weather_output <- "bc_weather_annual_metrics.csv"
merged_output <- "phenology_weather_merged.csv"
cor_output <- "phenology_weather_correlations.csv"

required_cols <- c("SPECIES", "SPNAME", "year_std", "julian_day", "Banded")
rpbo <- read.csv(rpbo_input, check.names = FALSE, stringsAsFactors = FALSE)
missing_cols <- setdiff(required_cols, names(rpbo))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in rpbo file:", paste(missing_cols, collapse = ", ")))
}

# Use only days with captures for phenology timing metrics.
rpbo <- rpbo[!is.na(rpbo$year_std) & !is.na(rpbo$julian_day) & !is.na(rpbo$Banded) & rpbo$Banded > 0, ]
if (nrow(rpbo) == 0) {
  stop("No rows with Banded > 0 were found in rpbo data.")
}

weighted_quantile <- function(x, w, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  out <- sapply(probs, function(p) x[which(cw >= p)[1]])
  names(out) <- paste0("q", probs * 100)
  out
}

weighted_mean_day <- function(x, w) {
  sum(x * w) / sum(w)
}

# 1) Extract focal species and split into arrival/departure groups.
species_split <- split(rpbo, rpbo$SPECIES)
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
if (n_species < 2) {
  stop("Need at least two species to define arrival vs departure groups.")
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

# 2) Compute per-species per-year phenology metrics from julian day.
key <- paste(rpbo$SPECIES, rpbo$year_std, sep = "__")
grouped <- split(rpbo, key)

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
    weighted_mean_julian_day = weighted_mean_day(d$julian_day, d$Banded),
    total_banded = sum(d$Banded),
    sampled_days = nrow(d),
    stringsAsFactors = FALSE
  )
})
phenology <- do.call(rbind, phenology_rows)
phenology <- phenology[order(phenology$SPECIES, phenology$year), ]
write.csv(phenology, phenology_output, row.names = FALSE, na = "")

# 3) Build annual weather metrics for each BC location file.
weather_files <- list.files(pattern = "_temp_clean\\.csv$")
if (length(weather_files) == 0) {
  stop("No *_temp_clean.csv files were found.")
}

weather_summaries <- lapply(weather_files, function(f) {
  station <- sub("_temp_clean\\.csv$", "", f)
  w <- read.csv(f, stringsAsFactors = FALSE)

  needed <- c("year", "julian_day", "min_temperature")
  missing_w <- setdiff(needed, names(w))
  if (length(missing_w) > 0) {
    stop(paste("Missing columns in", f, ":", paste(missing_w, collapse = ", ")))
  }

  years_of_interest <- sort(unique(phenology$year))
  w <- w[!is.na(w$year) & w$year %in% years_of_interest, ]

  by_year <- split(w, w$year)
  rows <- lapply(by_year, function(d) {
    data.frame(
      station = station,
      year = as.integer(d$year[1]),
      n_days = nrow(d),
      annual_mean_min_temp = mean(d$min_temperature, na.rm = TRUE),
      spring_mean_min_temp = mean(d$min_temperature[d$julian_day >= 60 & d$julian_day <= 151], na.rm = TRUE),
      summer_mean_min_temp = mean(d$min_temperature[d$julian_day >= 152 & d$julian_day <= 243], na.rm = TRUE),
      fall_mean_min_temp = mean(d$min_temperature[d$julian_day >= 244 & d$julian_day <= 304], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})

weather_annual <- do.call(rbind, weather_summaries)
weather_annual <- weather_annual[order(weather_annual$station, weather_annual$year), ]
write.csv(weather_annual, weather_output, row.names = FALSE, na = "")

# 4) Merge phenology + weather and calculate station-specific correlations.
merged <- merge(phenology, weather_annual, by = "year", all.x = TRUE, all.y = FALSE)
merged <- merged[order(merged$SPECIES, merged$station, merged$year), ]
write.csv(merged, merged_output, row.names = FALSE, na = "")

weather_metrics <- c(
  "annual_mean_min_temp",
  "spring_mean_min_temp",
  "summer_mean_min_temp",
  "fall_mean_min_temp"
)

combos <- split(merged, paste(merged$SPECIES, merged$station, sep = "__"))
cor_rows <- list()

for (nm in names(combos)) {
  d <- combos[[nm]]
  species <- d$SPECIES[1]
  spname <- d$SPNAME[1]
  station <- d$station[1]
  phase <- d$phase_group[1]
  phase_metric <- d$phase_metric[1]

  for (wm in weather_metrics) {
    keep <- !is.na(d$phase_julian_day) & !is.na(d[[wm]])
    dd <- d[keep, ]
    n <- nrow(dd)
    if (n < 8 || stats::sd(dd$phase_julian_day) == 0 || stats::sd(dd[[wm]]) == 0) {
      cor_rows[[length(cor_rows) + 1]] <- data.frame(
        SPECIES = species,
        SPNAME = spname,
        phase_group = phase,
        phase_metric = phase_metric,
        station = station,
        weather_metric = wm,
        n_years = n,
        cor_r = NA_real_,
        cor_p_value = NA_real_,
        slope_days_per_degC = NA_real_,
        lm_p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    ct <- suppressWarnings(cor.test(dd$phase_julian_day, dd[[wm]], method = "pearson"))
    fit <- lm(phase_julian_day ~ dd[[wm]], data = dd)
    slope <- unname(coef(fit)[2])
    fit_p <- summary(fit)$coefficients[2, 4]

    cor_rows[[length(cor_rows) + 1]] <- data.frame(
      SPECIES = species,
      SPNAME = spname,
      phase_group = phase,
      phase_metric = phase_metric,
      station = station,
      weather_metric = wm,
      n_years = n,
      cor_r = unname(ct$estimate),
      cor_p_value = ct$p.value,
      slope_days_per_degC = slope,
      lm_p_value = fit_p,
      stringsAsFactors = FALSE
    )
  }
}

cor_table <- do.call(rbind, cor_rows)
cor_table <- cor_table[order(-abs(cor_table$cor_r), cor_table$cor_p_value), ]
write.csv(cor_table, cor_output, row.names = FALSE, na = "")

cat("Species groups written:", species_group_output, "\n")
cat("Phenology metrics written:", phenology_output, "\n")
cat("Weather metrics written:", weather_output, "\n")
cat("Merged table written:", merged_output, "\n")
cat("Correlations written:", cor_output, "\n")
cat("Species:", n_species, "Years:", length(unique(phenology$year)), "\n")
cat("Rows in phenology:", nrow(phenology), "\n")
cat("Rows in weather annual:", nrow(weather_annual), "\n")
cat("Rows in correlation table:", nrow(cor_table), "\n")
