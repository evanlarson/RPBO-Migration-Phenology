options(stringsAsFactors = FALSE)

input_path <- "rpbo_clean.csv"
output_path <- "rpbo_clean_standardized.csv"
issues_path <- "rpbo_date_parse_issues.csv"

df <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!("DATE" %in% names(df))) {
  stop("Expected a DATE column in rpbo_clean.csv")
}

if (!("Year" %in% names(df))) {
  stop("Expected a Year column in rpbo_clean.csv")
}

normalize_year <- function(y_str) {
  y <- suppressWarnings(as.integer(y_str))
  if (is.na(y)) {
    return(NA_integer_)
  }
  if (nchar(y_str) == 2) {
    if (y <= 50) {
      return(2000L + y)
    }
    return(1900L + y)
  }
  y
}

safe_iso_date <- function(year, month, day) {
  if (any(is.na(c(year, month, day)))) {
    return(NA_character_)
  }
  if (month < 1 || month > 12 || day < 1 || day > 31) {
    return(NA_character_)
  }
  iso <- sprintf("%04d-%02d-%02d", year, month, day)
  parsed <- as.Date(iso, format = "%Y-%m-%d")
  if (is.na(parsed)) {
    return(NA_character_)
  }
  as.character(parsed)
}

date_raw <- trimws(df$DATE)
year_col <- suppressWarnings(as.integer(df$Year))

tokenize_slash <- function(x) {
  p <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(p) != 3) {
    return(c(NA_integer_, NA_integer_, NA_integer_))
  }
  c(
    suppressWarnings(as.integer(p[1])),
    suppressWarnings(as.integer(p[2])),
    normalize_year(p[3])
  )
}

is_slash <- grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}$", date_raw)

# Infer date-order convention for each Year using unambiguous slash dates only.
year_style <- setNames(rep(NA_character_, length(unique(year_col))), unique(year_col))
for (y in names(year_style)) {
  idx <- which(year_col == as.integer(y) & is_slash)
  if (!length(idx)) next

  parsed_tokens <- lapply(date_raw[idx], tokenize_slash)
  a <- vapply(parsed_tokens, `[`, integer(1), 1)
  b <- vapply(parsed_tokens, `[`, integer(1), 2)

  mdy_unambig <- sum(a >= 1 & a <= 12 & b > 12, na.rm = TRUE)
  dmy_unambig <- sum(a > 12 & b >= 1 & b <= 12, na.rm = TRUE)

  if (mdy_unambig > dmy_unambig) {
    year_style[[y]] <- "mdy"
  } else if (dmy_unambig > mdy_unambig) {
    year_style[[y]] <- "dmy"
  }
}

# Global fallback if a year has only ambiguous dates.
known_styles <- year_style[!is.na(year_style)]
global_style <- if (length(known_styles)) names(sort(table(known_styles), decreasing = TRUE))[1] else "mdy"

parse_result <- lapply(seq_along(date_raw), function(i) {
  raw <- date_raw[i]
  y_expected <- year_col[i]

  if (is.na(raw) || raw == "") {
    return(list(iso = NA_character_, note = "blank"))
  }

  # Numeric slash format (mixed MDY/DMY and 2/4-digit years)
  if (grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}$", raw)) {
    toks <- tokenize_slash(raw)
    a <- toks[1]
    b <- toks[2]
    y <- toks[3]

    if (any(is.na(c(a, b, y)))) {
      return(list(iso = NA_character_, note = "slash_non_numeric"))
    }

    if (a > 12 && b <= 12) {
      iso <- safe_iso_date(y, b, a)
      return(list(iso = iso, note = "slash_dmy_unambiguous"))
    }

    if (b > 12 && a <= 12) {
      iso <- safe_iso_date(y, a, b)
      return(list(iso = iso, note = "slash_mdy_unambiguous"))
    }

    if (a >= 1 && a <= 12 && b >= 1 && b <= 12) {
      style <- year_style[[as.character(y_expected)]]
      if (is.na(style)) {
        style <- global_style
      }

      if (style == "dmy") {
        iso <- safe_iso_date(y, b, a)
        return(list(iso = iso, note = "slash_ambiguous_inferred_dmy"))
      }

      iso <- safe_iso_date(y, a, b)
      return(list(iso = iso, note = "slash_ambiguous_inferred_mdy"))
    }

    return(list(iso = NA_character_, note = "slash_out_of_range"))
  }

  # "24-Sep-09" style
  if (grepl("^[0-9]{1,2}-[A-Za-z]{3}-[0-9]{2,4}$", raw)) {
    parsed <- as.Date(raw, format = "%d-%b-%y")
    if (is.na(parsed)) {
      parsed <- as.Date(raw, format = "%d-%b-%Y")
    }
    return(list(iso = as.character(parsed), note = "day_mon_text"))
  }

  # ISO style fallback
  if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", raw)) {
    parsed <- as.Date(raw, format = "%Y-%m-%d")
    return(list(iso = as.character(parsed), note = "iso"))
  }

  list(iso = NA_character_, note = "unknown_format")
})

date_iso <- vapply(parse_result, `[[`, character(1), "iso")
parse_note <- vapply(parse_result, `[[`, character(1), "note")

parsed_date <- as.Date(date_iso, format = "%Y-%m-%d")
year_std <- as.integer(format(parsed_date, "%Y"))
julian_day <- as.integer(format(parsed_date, "%j"))

# Cross-check: parsed year should match Year column where both exist.
year_mismatch <- !is.na(year_col) & !is.na(year_std) & (year_col != year_std)
if (any(year_mismatch)) {
  parse_note[year_mismatch] <- paste0(parse_note[year_mismatch], ";year_mismatch")
}

df$date_iso <- date_iso
df$year_std <- year_std
df$julian_day <- julian_day
df$date_parse_note <- parse_note

write.csv(df, output_path, row.names = FALSE, na = "")

issues <- df[is.na(df$date_iso) | grepl("year_mismatch", df$date_parse_note, fixed = TRUE), ]
write.csv(issues, issues_path, row.names = FALSE, na = "")

cat("Input rows:", nrow(df), "\n")
cat("Parsed dates:", sum(!is.na(df$date_iso)), "\n")
cat("Unparsed dates:", sum(is.na(df$date_iso)), "\n")
cat("Year mismatches:", sum(year_mismatch), "\n")
cat("Wrote:", output_path, "\n")
cat("Wrote issues:", issues_path, "\n")
