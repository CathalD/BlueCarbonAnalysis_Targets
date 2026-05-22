mod_guide_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "d-flex justify-content-between align-items-center mb-3",
      div(
        h4(class = "mb-0", "Run the analysis"),
        p(class = "text-muted mb-0 small",
          "Paste these commands into the RStudio Console in order. ",
          "Status updates when you click Refresh.")
      ),
      actionButton(ns("refresh"), "↺ Refresh status",
        class = "btn btn-outline-secondary btn-sm")
    ),
    uiOutput(ns("cards"))
  )
}

mod_guide_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {

    refresh_count <- reactiveVal(0L)
    observeEvent(input$refresh, refresh_count(refresh_count() + 1L))

    state <- reactive({
      refresh_count()   # re-read on manual refresh

      cfg <- read_project_config(project_root)

      config_ok <- file.exists(file.path(project_root, "blue_carbon_config.R"))
      locs_ok   <- file.exists(file.path(project_root,
        "Pre-Analysis Data Preparation", "data_raw", "core_locations.csv"))
      samp_ok   <- file.exists(file.path(project_root,
        "Pre-Analysis Data Preparation", "data_raw", "core_samples.csv"))

      covar_path <- cfg$COVARIATE_RASTER %||% ""
      has_raster <- nchar(trimws(covar_path)) > 0 && file.exists(covar_path)

      gee_project <- trimws(cfg$GEE_PROJECT %||% "")
      has_gee     <- nchar(gee_project) > 0

      list(
        setup_done  = config_ok && locs_ok && samp_ok,
        has_raster  = has_raster,
        has_gee     = has_gee,
        gee_project = gee_project,
        s1  = check_store_status(project_root, "_targets"),
        rf  = check_store_status(project_root, "_targets_rf"),
        pre = check_store_status(project_root, "_targets_preanalysis"),
        tl  = check_store_status(project_root, "_targets_transfer"),
        emb = check_store_status(project_root, "_targets_embedding")
      )
    })

    output$cards <- renderUI({
      s <- state()

      tagList(

        # ── Setup check ──────────────────────────────────────────────────────
        guide_setup_card(s),

        # ── Pipeline 1: Non-spatial ───────────────────────────────────────────
        pipeline_card(
          number      = 1,
          title       = "Non-spatial analysis",
          runtime     = "5–15 min",
          status      = s$s1,
          enabled     = s$setup_done,
          prereq_text = if (!s$setup_done) "Complete the Setup tab first.",
          description = paste0(
            "Harmonizes field cores to VM0033 standard depths, computes ",
            "per-stratum carbon stocks, generates exploratory plots. ",
            "Produces: ", code_span("reports/step1_nonspatial.html")),
          code = "targets::tar_make()"
        ),

        # ── Pipeline 2: RF spatial maps ───────────────────────────────────────
        pipeline_card(
          number      = 2,
          title       = "RF spatial prediction maps",
          runtime     = "5–10 min",
          status      = s$rf,
          enabled     = s$setup_done && s$has_raster,
          prereq_text = if (!s$has_raster)
            "Configure a covariate raster in the Setup tab (Step 3).",
          description = paste0(
            "Fits a random forest model using your field cores and the satellite ",
            "covariate raster. Produces a 25-m carbon stock map and ",
            code_span("reports/step3_random_forest.html")),
          code = 'targets::tar_make(\n  script = "_targets_rf.R",\n  store  = "_targets_rf"\n)'
        ),

        # ── Pipeline 3: GEE pre-analysis ─────────────────────────────────────
        pipeline_card(
          number      = 3,
          title       = "GEE global covariate extraction",
          runtime     = "~60 min",
          status      = s$pre,
          enabled     = s$setup_done && s$has_gee,
          prereq_text = if (!s$has_gee)
            "Enter a GEE project ID in the Setup tab (Step 1) to enable this pipeline.",
          description = paste0(
            "Extracts 26 satellite covariates at ~952 global wetland cores via ",
            "Google Earth Engine. Run once — required before Pipelines 4 and 5."),
          extra = if (s$setup_done && s$has_gee) {
            div(class = "alert alert-info py-2 px-3 mb-2",
              tags$strong("One-time GEE authentication — run this in RStudio first:"),
              tags$pre(class = "code-block mb-0 mt-1",
                paste0("library(rgee)\n",
                       "ee_Initialize(user = \"your.email@gmail.com\", drive = TRUE)"))
            )
          },
          code = 'targets::tar_make(\n  script = "_targets_preanalysis.R",\n  store  = "_targets_preanalysis"\n)'
        ),

        # ── Pipeline 4: Wadoux TL ─────────────────────────────────────────────
        pipeline_card(
          number      = 4,
          title       = "Wadoux transfer learning",
          runtime     = "~15 min",
          status      = s$tl,
          enabled     = s$pre == "complete",
          prereq_text = if (s$pre != "complete")
            "Run Pipeline 3 (GEE pre-analysis) first.",
          description = paste0(
            "Weights ~952 global wetland cores by similarity to your site ",
            "and trains a bias-corrected RF with bootstrap uncertainty. ",
            "Produces: ", code_span("reports/step4_transfer_learning.html")),
          code = 'targets::tar_make(\n  script = "_targets_transfer.R",\n  store  = "_targets_transfer"\n)'
        ),

        # ── Pipeline 5: Embedding TL ──────────────────────────────────────────
        pipeline_card(
          number      = 5,
          title       = "Embedding transfer learning",
          runtime     = "~30 min",
          status      = s$emb,
          enabled     = s$pre == "complete",
          prereq_text = if (s$pre != "complete")
            "Run Pipeline 3 (GEE pre-analysis) first.",
          description = paste0(
            "Uses Google's 64-d satellite foundation model for site similarity ",
            "instead of a domain classifier. Produces maps comparable to Pipeline 4 ",
            "plus ", code_span("reports/step5_embedding_tl.html"), " with a comparison table."),
          code = 'targets::tar_make(\n  script = "_targets_embedding.R",\n  store  = "_targets_embedding"\n)'
        ),

        # ── Troubleshooting ───────────────────────────────────────────────────
        bslib::card(
          class = "mt-3 border-0 bg-light",
          bslib::card_body(
            tags$strong("If a pipeline errors:"),
            tags$pre(class = "code-block mt-1 mb-0",
'targets::tar_meta() |>
  dplyr::filter(!is.na(error)) |>
  dplyr::select(name, error)')
          )
        )
      )
    })
  })
}

# ── Helpers ────────────────────────────────────────────────────────────────────

check_store_status <- function(project_root, store_name) {
  store_path <- file.path(project_root, store_name)
  if (!dir.exists(store_path)) return("not_run")
  meta <- tryCatch(
    targets::tar_meta(store = store_path),
    error = function(e) NULL
  )
  if (is.null(meta) || nrow(meta) == 0) return("not_run")
  if (any(!is.na(meta$error))) return("error")
  "complete"
}

read_project_config <- function(project_root) {
  config_path <- file.path(project_root, "blue_carbon_config.R")
  if (!file.exists(config_path)) return(list())
  e <- new.env(parent = emptyenv())
  tryCatch(
    { source(config_path, local = e, echo = FALSE); as.list(e) },
    error = function(x) list()
  )
}

code_span <- function(text) tags$code(text)

guide_setup_card <- function(s) {
  items <- list(
    list(
      label = "Configuration file (blue_carbon_config.R)",
      ok    = file.exists(
        file.path(
          tryCatch(normalizePath(file.path(getwd(), "..")), error = function(e) "."),
          "blue_carbon_config.R"
        )
      )
    )
  )
  # Use s$setup_done as combined indicator
  cls <- if (s$setup_done) "border-success" else "border-warning"
  icon <- if (s$setup_done) span(class = "badge bg-success", "✓ Ready") else
    span(class = "badge bg-warning text-dark", "Incomplete")

  bslib::card(
    class = paste("mb-3", cls),
    style = "border-left-width: 4px;",
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center",
      tags$strong("Setup"),
      icon
    ),
    if (!s$setup_done) {
      bslib::card_body(
        p(class = "mb-0 text-muted",
          "Complete the ", tags$strong("Setup tab"), " and click ",
          tags$strong("Save Setup Files"), " before running any pipelines.")
      )
    }
  )
}

status_badge <- function(status) {
  switch(status,
    not_run  = span(class = "badge bg-secondary", "Not run"),
    complete = span(class = "badge bg-success",   "✓ Complete"),
    error    = span(class = "badge bg-danger",    "✗ Error — check tar_meta()"),
    span(class = "badge bg-secondary", "Not run")
  )
}

pipeline_card <- function(number, title, runtime, status, enabled,
                           prereq_text = NULL, description = NULL,
                           extra = NULL, code) {
  card_border <- switch(status,
    complete = "border-success",
    error    = "border-danger",
    ""
  )

  bslib::card(
    class = paste("mb-3", card_border, if (!enabled) "opacity-60"),
    style = if (status %in% c("complete", "error")) "border-left-width: 4px;" else NULL,
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center",
      div(class = "d-flex align-items-center gap-2",
        tags$strong(paste0("Pipeline ", number, " — ", title)),
        status_badge(status)
      ),
      span(class = "badge bg-light text-dark border", runtime)
    ),
    bslib::card_body(
      if (!is.null(prereq_text) && !enabled) {
        div(class = "alert alert-warning py-2 px-3 mb-2 small",
          tags$strong("Prerequisite: "), prereq_text)
      },
      if (!is.null(description)) p(class = "text-muted small mb-2", description),
      extra,
      tags$pre(class = "code-block mb-0", code)
    )
  )
}
