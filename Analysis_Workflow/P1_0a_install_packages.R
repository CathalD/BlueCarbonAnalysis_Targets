# ============================================================================
# MODULE 00a: BLUE CARBON - PACKAGE INSTALLATION (SIMPLIFIED)
# ============================================================================
# PURPOSE: Install required R packages with robust error handling
# USAGE: Run this FIRST, then run 00b_setup_directories.R
# ============================================================================

cat("\n========================================\n")
cat("BLUE CARBON - PACKAGE INSTALLATION\n")
cat("========================================\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set options for clean installation
options(
  repos = c(CRAN = "https://cloud.r-project.org/"),
  timeout = 300,  # 5 minute timeout per package
  install.packages.compile.from.source = "never"  # Binary only!
)

# ============================================================================
# CHECK R VERSION
# ============================================================================

cat("Checking R version...\n")
r_version <- paste(R.version$major, R.version$minor, sep = ".")
cat(sprintf("  R version: %s\n", r_version))

if (as.numeric(R.version$major) < 4) {
  cat("  ⚠️  R version 4.0+ recommended\n\n")
} else {
  cat("  ✓ R version OK\n\n")
}

# ============================================================================
# DEFINE REQUIRED PACKAGES
# ============================================================================

# Core packages (always required)
required_packages <- c(
  # Data manipulation
  "dplyr", "tidyr", "readr",

  # Visualization
  "ggplot2", "gridExtra",

  # Spatial analysis
  "sf", "terra", "gstat",

  # Statistical modeling
  "randomForest", "mgcv", "boot",

  # Blue carbon specific
  "CAST",

  # Reporting
  "openxlsx"
)

# Optional packages (nice to have)
optional_packages <- c(
  "raster",      # Legacy spatial (terra is newer)
  "spdep",       # Spatial dependence
  "corrplot",    # Correlation plots
  "lubridate",   # Date handling
  "aqp",         # Soil profile analysis
  "ithir",       # TRUE equal-area spline for depth harmonization
  "caret",       # ML framework
  "knitr",       # Report generation
  "viridis",     # Color palettes
  "isotree"      # Isolation forest for outlier detection
)

cat(sprintf("Required packages: %d\n", length(required_packages)))
cat(sprintf("Optional packages: %d\n\n", length(optional_packages)))

# ============================================================================
# SIMPLE INSTALLATION FUNCTION
# ============================================================================

install_pkg <- function(pkg, verbose = TRUE) {

  # Check if already installed
  if (requireNamespace(pkg, quietly = TRUE)) {
    if (verbose) cat(sprintf("  ✓ %s (already installed)\n", pkg))
    return(TRUE)
  }

  # Try to install
  if (verbose) cat(sprintf("  Installing %s... ", pkg))

  success <- tryCatch({
    # Suppress all output during installation
    suppressWarnings(
      suppressMessages(
        install.packages(pkg,
                        dependencies = TRUE,
                        quiet = TRUE,
                        type = "binary",  # Binary only - no compilation!
                        repos = "https://cloud.r-project.org/")
      )
    )

    # Verify installation
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (verbose) cat("✓\n")
      TRUE
    } else {
      if (verbose) cat("✗ (verification failed)\n")
      FALSE
    }

  }, error = function(e) {
    if (verbose) cat(sprintf("✗ (error: %s)\n", e$message))
    FALSE
  }, warning = function(w) {
    # Treat warnings as success if package loads
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (verbose) cat("✓\n")
      TRUE
    } else {
      if (verbose) cat("✗\n")
      FALSE
    }
  })

  return(success)
}

# ============================================================================
# INSTALL REQUIRED PACKAGES
# ============================================================================

cat("\n========================================\n")
cat("INSTALLING REQUIRED PACKAGES\n")
cat("========================================\n\n")

cat("Note: Using binary packages only (no compilation)\n")
cat("This is faster and more reliable.\n\n")

required_results <- logical(length(required_packages))
names(required_results) <- required_packages

for (i in seq_along(required_packages)) {
  pkg <- required_packages[i]
  cat(sprintf("[%d/%d] ", i, length(required_packages)))
  required_results[pkg] <- install_pkg(pkg, verbose = TRUE)
  Sys.sleep(0.1)  # Small pause between packages
}

required_success <- sum(required_results)
required_total <- length(required_packages)

cat(sprintf("\nRequired packages: %d/%d installed (%.1f%%)\n\n",
            required_success, required_total,
            100 * required_success / required_total))

# ============================================================================
# INSTALL OPTIONAL PACKAGES
# ============================================================================

cat("\n========================================\n")
cat("INSTALLING OPTIONAL PACKAGES\n")
cat("========================================\n\n")

cat("Optional packages enhance functionality but aren't critical.\n\n")

optional_results <- logical(length(optional_packages))
names(optional_results) <- optional_packages

for (i in seq_along(optional_packages)) {
  pkg <- optional_packages[i]
  cat(sprintf("[%d/%d] ", i, length(optional_packages)))
  optional_results[pkg] <- install_pkg(pkg, verbose = TRUE)
  Sys.sleep(0.1)
}

optional_success <- sum(optional_results)
optional_total <- length(optional_packages)

cat(sprintf("\nOptional packages: %d/%d installed (%.1f%%)\n\n",
            optional_success, optional_total,
            100 * optional_success / optional_total))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("INSTALLATION SUMMARY\n")
cat("========================================\n\n")

# Calculate overall success
total_attempted <- required_total + optional_total
total_success <- required_success + optional_success
overall_rate <- 100 * total_success / total_attempted

cat(sprintf("Total packages: %d/%d installed (%.1f%%)\n\n",
            total_success, total_attempted, overall_rate))

# List missing required packages
missing_required <- required_packages[!required_results]
if (length(missing_required) > 0) {
  cat("⚠️  MISSING REQUIRED PACKAGES:\n")
  for (pkg in missing_required) {
    cat(sprintf("  ✗ %s\n", pkg))
  }
  cat("\n")
}

# List missing optional packages
missing_optional <- optional_packages[!optional_results]
if (length(missing_optional) > 0) {
  cat("Missing optional packages (not critical):\n")
  for (pkg in missing_optional) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\n")
}

# Check critical packages for workflow
critical_packages <- c("dplyr", "ggplot2", "sf", "terra", "randomForest", "openxlsx")
critical_status <- sapply(critical_packages,
                         function(p) requireNamespace(p, quietly = TRUE))

cat("Critical packages for workflow:\n")
for (i in seq_along(critical_packages)) {
  pkg <- critical_packages[i]
  status <- if (critical_status[i]) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, pkg))
}
cat("\n")

# ============================================================================
# NEXT STEPS
# ============================================================================

cat("========================================\n")

if (required_success == required_total) {
  cat("✓✓✓ SUCCESS! All required packages installed.\n\n")
  cat("Next step:\n")
  cat("  source('00b_setup_directories.R')\n\n")

} else if (required_success >= 0.9 * required_total) {
  cat("✓ MOSTLY COMPLETE - Workflow should work.\n\n")

  if (length(missing_required) > 0) {
    cat("Try installing missing packages manually:\n")
    for (pkg in missing_required) {
      cat(sprintf("  install.packages('%s', type = 'binary')\n", pkg))
    }
    cat("\n")
  }

  cat("Next step:\n")
  cat("  source('00b_setup_directories.R')\n\n")

} else {
  cat("⚠️  INCOMPLETE INSTALLATION\n\n")

  cat("Manual installation required for missing packages.\n")
  cat("Try these commands one at a time:\n\n")

  for (pkg in missing_required) {
    cat(sprintf("install.packages('%s', type = 'binary')\n", pkg))
  }
  cat("\n")

  cat("For spatial packages, you may need system libraries:\n")
  cat("  Mac:    brew install gdal proj geos\n")
  cat("  Ubuntu: sudo apt-get install gdal-bin libgdal-dev libproj-dev\n\n")
}

cat("Done! 🌊\n\n")

# ============================================================================
# SAVE INSTALLATION RECORD
# ============================================================================

if (!dir.exists("logs")) {
  dir.create("logs", recursive = TRUE, showWarnings = FALSE)
}

install_record <- data.frame(
  package = c(required_packages, optional_packages),
  category = c(rep("required", length(required_packages)),
               rep("optional", length(optional_packages))),
  installed = c(required_results, optional_results),
  date = Sys.Date(),
  stringsAsFactors = FALSE
)

write.csv(install_record,
          file.path("logs", paste0("package_install_", Sys.Date(), ".csv")),
          row.names = FALSE)

cat("Installation record saved to: logs/package_install_", Sys.Date(), ".csv\n\n", sep = "")
