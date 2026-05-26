# .Rprofile — loaded automatically when opening this RStudio project
library(targets)
library(tarchetypes)

# Main pipeline (non-spatial: Steps 1–2 + report)
tm   <- function(...) targets::tar_make(...)
tv   <- function()    targets::tar_visnetwork()
tl   <- function(x)  targets::tar_load(!!rlang::ensym(x))
tr   <- function(x)  targets::tar_read(x)

# RF pipeline (Step 3: spatial prediction maps)
tmrf <- function()   targets::tar_make(script = "_targets_rf.R",         store = "_targets_rf")

# BlueCarbon compaction correction + stock comparison
tmbc <- function()   targets::tar_make(script = "_targets_bluecarbon.R", store = "_targets_bluecarbon")

# Sequestration rates (requires core_chronology.csv)
tmsr <- function()   targets::tar_make(script = "_targets_seqrates.R",   store = "_targets_seqrates")

app  <- function()   shiny::runApp("shiny", launch.browser = TRUE)

message("Shortcuts: tm() | tmrf() | tmbc() | tmsr() | app()")
