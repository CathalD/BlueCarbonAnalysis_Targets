validate_locations_csv <- function(df, valid_strata = NULL) {
  errors   <- character(0)
  warnings <- character(0)

  required_cols <- c("core_id", "latitude", "longitude", "stratum")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing, collapse = ", ")))
    return(list(errors = errors, warnings = warnings))
  }

  if (nrow(df) == 0) {
    errors <- c(errors, "File is empty — no rows found.")
    return(list(errors = errors, warnings = warnings))
  }

  if (any(is.na(df$core_id) | trimws(as.character(df$core_id)) == "")) {
    errors <- c(errors, "core_id contains empty or missing values.")
  }

  dups <- df$core_id[duplicated(df$core_id)]
  if (length(dups) > 0) {
    errors <- c(errors, paste("Duplicate core_id values:", paste(unique(dups), collapse = ", ")))
  }

  lat <- suppressWarnings(as.numeric(df$latitude))
  if (anyNA(lat)) {
    errors <- c(errors, "latitude contains non-numeric values.")
  } else if (any(lat < -90 | lat > 90)) {
    errors <- c(errors, "latitude values must be between -90 and 90 (decimal degrees).")
  }

  lon <- suppressWarnings(as.numeric(df$longitude))
  if (anyNA(lon)) {
    errors <- c(errors, "longitude contains non-numeric values.")
  } else if (any(lon < -180 | lon > 180)) {
    errors <- c(errors, "longitude values must be between -180 and 180 (decimal degrees).")
  }

  if (!is.null(valid_strata) && length(valid_strata) > 0) {
    unknown <- setdiff(unique(as.character(df$stratum)), valid_strata)
    if (length(unknown) > 0) {
      warnings <- c(warnings, paste0(
        "Stratum value(s) not in VALID_STRATA (",
        paste(valid_strata, collapse = ", "), "): ",
        paste(unknown, collapse = ", "),
        ". These cores will be excluded from analysis."
      ))
    }
  }

  list(errors = errors, warnings = warnings)
}

validate_samples_csv <- function(df) {
  errors   <- character(0)
  warnings <- character(0)

  required_cols <- c("core_id", "depth_top_cm", "depth_bottom_cm", "soc_g_kg")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing, collapse = ", ")))
    return(list(errors = errors, warnings = warnings))
  }

  if (nrow(df) == 0) {
    errors <- c(errors, "File is empty — no rows found.")
    return(list(errors = errors, warnings = warnings))
  }

  for (col in c("depth_top_cm", "depth_bottom_cm", "soc_g_kg")) {
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (anyNA(vals)) {
      errors <- c(errors, paste0("'", col, "' contains non-numeric values."))
    }
  }

  top <- suppressWarnings(as.numeric(df$depth_top_cm))
  bot <- suppressWarnings(as.numeric(df$depth_bottom_cm))
  if (!anyNA(top) && !anyNA(bot)) {
    if (any(top >= bot)) {
      errors <- c(errors, "Some rows have depth_top_cm ≥ depth_bottom_cm (top must be less than bottom).")
    }
    if (any(top < 0)) {
      errors <- c(errors, "depth_top_cm contains negative values.")
    }
  }

  soc <- suppressWarnings(as.numeric(df$soc_g_kg))
  if (!anyNA(soc) && any(soc < 0 | soc > 500)) {
    warnings <- c(warnings, paste0(
      "Some SOC values are outside 0–500 g/kg. ",
      "Check for unit errors — SOC should be in g C per kg dry soil."
    ))
  }

  if ("bulk_density_g_cm3" %in% names(df)) {
    bd <- suppressWarnings(as.numeric(df$bulk_density_g_cm3))
    n_missing <- sum(is.na(bd))
    if (n_missing > 0) {
      warnings <- c(warnings, paste0(
        n_missing, " sample(s) are missing bulk_density_g_cm3 — ",
        "the ecosystem default (0.8 g/cm³) will be used for those rows."
      ))
    }
    bd_ok <- bd[!is.na(bd)]
    if (length(bd_ok) > 0 && any(bd_ok < 0.1 | bd_ok > 3.0)) {
      warnings <- c(warnings, "Some bulk density values are outside the normal range (0.1–3.0 g/cm³).")
    }
  } else {
    warnings <- c(warnings, paste0(
      "Column 'bulk_density_g_cm3' not found — ",
      "the ecosystem default (0.8 g/cm³) will be used for all samples."
    ))
  }

  list(errors = errors, warnings = warnings)
}
