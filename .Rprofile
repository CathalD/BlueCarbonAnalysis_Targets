# .Rprofile — loaded automatically when opening this RStudio project
library(targets)
library(tarchetypes)

tm  <- function(...) targets::tar_make(...)
tv  <- function()    targets::tar_visnetwork()
tl  <- function(x)  targets::tar_load(!!rlang::ensym(x))
tr  <- function(x)  targets::tar_read(x)
tm1 <- function()   targets::tar_make(names = c(
  "locations_file", "samples_file", "covar_file", "cfg",
  "cores_raw", "eda_plots", "cores_harmonized",
  "stratum_summary", "step2_extrapolation",
  "rf_data", "rf_models", "rf_rasters", "rf_importance_plot", "rf_maps",
  "report_nonspatial", "report_rf"
))

message("targets loaded. Use tm(), tv(), tl(target_name), tm1() for full pipeline.")
