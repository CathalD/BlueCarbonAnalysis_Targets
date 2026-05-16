# .Rprofile — loaded automatically when opening this RStudio project
library(targets)
library(tarchetypes)

# Short aliases — type these in the R console during development
tm  <- function(...) targets::tar_make(...)
tv  <- function()    targets::tar_visnetwork()
tl  <- function(x)  targets::tar_load(!!rlang::ensym(x))
tr  <- function(x)  targets::tar_read(x)
tm1 <- function()   targets::tar_make(names = c(
  "locations_file", "samples_file", "cfg",
  "cores_raw", "cores_clean", "eda_plots",
  "cores_harmonized", "stratum_summary", "report_step1"
))

message("targets loaded. Use tm(), tv(), tl(target_name), tm1() for Step 1 only.")
