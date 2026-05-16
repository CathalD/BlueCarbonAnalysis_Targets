# R/config.R
# ============================================================================
# PURPOSE: Wrap blue_carbon_config.R into a function that returns a list.
#
# WHY A LIST?
#   In targets, each target is a self-contained function call. Functions cannot
#   rely on objects in .GlobalEnv the way source() scripts can. By returning a
#   named list, every downstream function receives only what it needs explicitly.
#
# HOW targets USES THIS:
#   In _targets.R:
#     tar_target(cfg, load_config())
#   targets calls load_config() once, stores the result as "cfg", and passes
#   it to any downstream target that lists cfg as an argument.
#
# ACCESSING VALUES:
#   cfg$PROJECT_NAME
#   cfg$VM0033_DEPTH_MIDPOINTS        # c(7.5, 22.5, 40, 75)
#   cfg$BD_DEFAULTS                   # list(IM = 0.8, NM = 0.8, MF = 0.8)
#   cfg$QC_SOC_MIN / cfg$QC_SOC_MAX   # 0 / 500
#   cfg$QC_BD_MIN  / cfg$QC_BD_MAX    # 0.1 / 3.0
# ============================================================================
load_config <- function(config_path = "blue_carbon_config.R") {
  # Source the config into a fresh local environment.
  # This does NOT pollute .GlobalEnv — all assignments stay inside env.
  env <- new.env(parent = baseenv())
  source(config_path, local = env)
  # Return everything as a named list.
  as.list(env)
}
