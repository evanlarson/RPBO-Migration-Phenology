options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

in_results <- "model_reduced_hypothesis_results.csv"
in_phen <- "rpbo_focal_species_phenology_jul_oct.csv"
in_weather <- "bc_weather_region_monthly_wide_jul_oct.csv"
in_net_hours <- "nethours.xlsx"
out_file <- "model_reduced_significant_loo_sensitivity.csv"

res <- read.csv(in_results, stringsAsFactors = FALSE)
phen <- read.csv(in_phen, stringsAsFactors = FALSE)
weather <- read.csv(in_weather, stringsAsFactors = FALSE)

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
  net_col <- names(d)[grepl("net", nms) & grepl("hour", nms)][1]
  if (is.na(year_col) || is.na(net_col)) {
    stop("Net-hours file must include year and net-hours columns.")
  }
  data.frame(
    year = suppressWarnings(as.integer(d[[year_col]])),
    net_hours = suppressWarnings(as.numeric(d[[net_col]])),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(year), !is.na(net_hours), net_hours > 0) %>%
    group_by(year) %>%
    summarize(net_hours = mean(net_hours, na.rm = TRUE), .groups = "drop") %>%
    mutate(net_hours_100 = net_hours / 100)
}

net_hours <- load_net_hours(in_net_hours)

metrics_lookup <- data.frame(
  metric_code = c("first_julian_day", "q10_julian_day", "q50_julian_day", "q90_julian_day", "last_julian_day"),
  metric_label = c("first_detection", "p10", "p50", "p90", "last_detection"),
  stringsAsFactors = FALSE
)

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

sig <- res %>%
  filter(
    temp_p_fdr_global_glm < 0.10 |
      temp_p_fdr_global_permutation < 0.10 |
      temp_p_fdr_global_gam < 0.10
  ) %>%
  filter(
    (phase_group == "departure" & metric_label != "first_detection") |
      (phase_group == "fall_arrival" & metric_label != "last_detection")
  )

if (nrow(sig) == 0) {
  write.csv(data.frame(), out_file, row.names = FALSE)
  cat("No FDR-significant rows found; wrote empty file:", out_file, "\n")
  quit(save = "no", status = 0)
}

rows <- vector("list", nrow(sig))
for (i in seq_len(nrow(sig))) {
  s <- sig[i, ]
  temp_var <- paste0(s$temp_month, "_mean_min_temp")

  d_sp <- phen_long %>%
    filter(SPECIES == s$SPECIES, metric_label == s$metric_label) %>%
    select(year, timing_julian_day)
  d_rg <- weather %>%
    filter(region == s$region) %>%
    select(year, all_of(temp_var))

  dd <- inner_join(d_sp, d_rg, by = "year") %>%
    rename(temp = all_of(temp_var)) %>%
    left_join(net_hours %>% select(year, net_hours_100), by = "year") %>%
    filter(!is.na(timing_julian_day), !is.na(temp), !is.na(net_hours_100), !is.na(year)) %>%
    arrange(year)

  fit <- lm(timing_julian_day ~ year + net_hours_100 + temp, data = dd)
  beta <- unname(coef(fit)["temp"])
  p <- summary(fit)$coefficients["temp", "Pr(>|t|)"]

  loo_beta <- rep(NA_real_, nrow(dd))
  loo_p <- rep(NA_real_, nrow(dd))
  for (j in seq_len(nrow(dd))) {
    d2 <- dd[-j, ]
    f2 <- lm(timing_julian_day ~ year + net_hours_100 + temp, data = d2)
    loo_beta[j] <- unname(coef(f2)["temp"])
    loo_p[j] <- summary(f2)$coefficients["temp", "Pr(>|t|)"]
  }

  rows[[i]] <- data.frame(
    SPECIES = s$SPECIES,
    SPNAME = s$SPNAME,
    metric_label = s$metric_label,
    region = s$region,
    temp_month = s$temp_month,
    n_years = nrow(dd),
    beta_temp_full = beta,
    p_temp_full = p,
    sign_consistency = mean(sign(loo_beta) == sign(beta), na.rm = TRUE),
    loo_p_lt_0_05_rate = mean(loo_p < 0.05, na.rm = TRUE),
    loo_beta_min = min(loo_beta, na.rm = TRUE),
    loo_beta_max = max(loo_beta, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

out <- bind_rows(rows) %>% arrange(p_temp_full)
write.csv(out, out_file, row.names = FALSE, na = "")
cat("Wrote:", out_file, "\n")
