mod_finish_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Choose your analyses and save."),
      p("The non-spatial analysis always runs. Spatial and transfer learning options ",
        "are available when a raster and GEE project are configured.")
    ),

    # ── What you'll get ────────────────────────────────────────────────────
    uiOutput(ns("analysis_menu")),

    # ── Configuration summary (collapsible) ───────────────────────────────
    tags$details(class = "mb-3",
      tags$summary(class = "text-muted small",
        style = "cursor:pointer;", "View configuration summary"),
      uiOutput(ns("config_summary"))
    ),

    # ── Save button + result ───────────────────────────────────────────────
    uiOutput(ns("save_section"))
  )
}

mod_finish_server <- function(id, setup_state, data_state, raster_state, project_root) {
  moduleServer(id, function(input, output, session) {
    ns      <- session$ns
    saved   <- reactiveVal(FALSE)
    save_err <- reactiveVal(NULL)

    # ── Derived capability flags ───────────────────────────────────────────
    has_raster <- reactive({ isTRUE(raster_state()$has_raster) })
    has_gee    <- reactive({ nchar(trimws(setup_state()$gee_project)) > 0 })

    # ── Analysis menu ──────────────────────────────────────────────────────
    output$analysis_menu <- renderUI({
      r <- has_raster()
      g <- has_gee()

      bslib::card(
        bslib::card_header("What will be produced"),
        bslib::card_body(
          pipeline_row(
            ns       = ns,
            input_id = NULL,          # always runs, no checkbox
            checked  = TRUE,
            enabled  = TRUE,
            label    = "Non-spatial analysis",
            badge    = "Always included",
            badge_cls = "bg-success",
            items    = c(
              "Depth harmonization to VM0033 standard intervals",
              "Per-stratum carbon stock table (kg C/m²) — core VM0033 reporting values",
              "Exploratory plots: depth profiles, spatial map, stock distributions",
              "HTML report: reports/step1_nonspatial.html"
            )
          ),

          pipeline_row(
            ns        = ns,
            input_id  = "run_rf",
            checked   = r,
            enabled   = r,
            label     = "Local random forest — spatial prediction map",
            badge     = if (r) "Raster configured" else "Requires raster",
            badge_cls = if (r) "bg-primary" else "bg-secondary",
            items     = c(
              "25-m resolution carbon stock map across the whole site",
              "Variable importance: which satellite bands drive predictions",
              "Leave-one-core-out cross-validation accuracy (R², RMSE)",
              "HTML report: reports/step3_random_forest.html"
            )
          ),

          pipeline_row(
            ns        = ns,
            input_id  = "run_tl",
            checked   = FALSE,
            enabled   = r && g,
            label     = "Transfer learning — Wadoux method",
            badge     = if (r && g) "Optional" else "Requires raster + GEE",
            badge_cls = if (r && g) "bg-info text-dark" else "bg-secondary",
            items     = c(
              "Borrows signal from ~952 global wetland cores weighted by site similarity",
              "Bias-corrected predictions + 500-replicate bootstrap uncertainty",
              "HTML report: reports/step4_transfer_learning.html"
            )
          ),

          pipeline_row(
            ns        = ns,
            input_id  = "run_emb",
            checked   = FALSE,
            enabled   = r && g,
            label     = "Transfer learning — Embedding method",
            badge     = if (r && g) "Optional" else "Requires raster + GEE",
            badge_cls = if (r && g) "bg-info text-dark" else "bg-secondary",
            items     = c(
              "Uses Google's 64-d satellite foundation model instead of a classifier",
              "Produces same maps as Wadoux TL but with a different similarity measure",
              "HTML report: reports/step5_embedding_tl.html (includes Wadoux comparison)"
            )
          )
        )
      )
    })

    # ── Config summary (collapsed) ─────────────────────────────────────────
    output$config_summary <- renderUI({
      s <- setup_state()
      d <- data_state()
      r <- raster_state()
      n_cores   <- if (!is.null(d$locations)) nrow(d$locations) else 0
      n_samples <- if (!is.null(d$samples))   nrow(d$samples)   else 0

      div(class = "mt-2",
        tags$table(class = "table table-sm table-borderless mb-0",
          tags$tbody(
            summary_row("Project",      s$project_name),
            summary_row("Location",     s$project_location),
            summary_row("Year",         as.character(s$monitoring_year)),
            summary_row("Strata",       paste(s$valid_strata, collapse = ", ")),
            summary_row("GEE project",
              if (nchar(s$gee_project) > 0) s$gee_project
              else tags$em("Not set")),
            summary_row("Cores",        paste0(n_cores, " locations, ", n_samples, " depth intervals")),
            summary_row("Raster",
              if (r$has_raster) tags$span(class = "text-success", "✓ ", basename(r$raster_path))
              else tags$span(class = "text-muted", "Not configured")),
            summary_row("AOI boundary",
              if (!is.null(r$aoi_source)) tags$span(class = "text-success", "✓ Uploaded")
              else tags$span(class = "text-muted", "Not provided"))
          )
        )
      )
    })

    # ── Save section ───────────────────────────────────────────────────────
    output$save_section <- renderUI({
      if (saved()) return(commands_ui(ns, setup_state, raster_state, input))

      tagList(
        if (!is.null(save_err())) {
          div(class = "alert alert-danger", tags$strong("Save failed: "), save_err())
        },
        actionButton(ns("save_btn"), "Save Setup Files",
          class = "btn btn-success btn-lg",
          icon  = icon("save"))
      )
    })

    # ── Save handler ───────────────────────────────────────────────────────
    observeEvent(input$save_btn, {
      s <- setup_state()
      d <- data_state()
      r <- raster_state()

      result <- tryCatch({
        data_raw <- file.path(project_root,
          "Pre-Analysis Data Preparation", "data_raw")
        dir.create(data_raw, showWarnings = FALSE, recursive = TRUE)

        readr::write_csv(d$locations, file.path(data_raw, "core_locations.csv"))
        readr::write_csv(d$samples,   file.path(data_raw, "core_samples.csv"))

        if (!is.null(r$aoi_source) && !is.null(r$aoi_dest)) {
          dir.create(dirname(r$aoi_dest), showWarnings = FALSE, recursive = TRUE)
          file.copy(r$aoi_source, r$aoi_dest, overwrite = TRUE)
        }

        write_config(
          project_name     = s$project_name,
          project_location = s$project_location,
          monitoring_year  = s$monitoring_year,
          valid_strata     = s$valid_strata,
          gee_project      = s$gee_project,
          covariate_raster = r$raster_path,
          aoi_file         = r$aoi_dest,
          output_path      = file.path(project_root, "blue_carbon_config.R")
        )
        TRUE
      }, error = function(e) { save_err(conditionMessage(e)); FALSE })

      if (result) { saved(TRUE); save_err(NULL) }
    })
  })
}

# ── Commands panel (shown after save) ─────────────────────────────────────
commands_ui <- function(ns, setup_state, raster_state, input) {
  s <- setup_state()
  r <- raster_state()
  has_r <- isTRUE(r$has_raster)
  has_g <- nchar(trimws(s$gee_project)) > 0

  run_rf  <- isTRUE(input$run_rf)  && has_r
  run_tl  <- isTRUE(input$run_tl)  && has_r && has_g
  run_emb <- isTRUE(input$run_emb) && has_r && has_g

  tagList(
    div(class = "alert alert-success mt-3",
      tags$strong("✓ Setup files saved."),
      " Close this app, open the project in RStudio, and run the commands below."
    ),

    bslib::card(
      class = "mt-3 border-primary",
      bslib::card_header(class = "bg-primary text-white",
        "Run the analysis — paste into the RStudio Console"),
      bslib::card_body(

        # Always: Pipeline 1
        div(class = "pipeline-block",
          tags$p(tags$strong("Step 1 — Core data processing"),
            tags$span(class = "badge bg-success ms-2", "always runs")),
          if (!has_r) {
            p(class = "text-muted small",
              "Produces: depth harmonization, stratum summary, VM0033 estimates, ",
              "and ", code("reports/step1_nonspatial.html"))
          } else {
            p(class = "text-muted small",
              "Produces: non-spatial report + RF spatial map")
          },
          tags$pre(class = "code-block", "targets::tar_make()")
        ),

        # GEE auth block (shown if any GEE pipeline selected)
        if ((run_tl || run_emb) && has_g) {
          tagList(
            hr(),
            div(class = "alert alert-info py-2",
              tags$strong("GEE authentication required for steps 2–4."),
              " Run once — a browser will open to log in:"
            ),
            tags$pre(class = "code-block",
              paste0('library(rgee)\n',
                     'ee_Initialize(user = "your.email@gmail.com", drive = TRUE)'))
          )
        },

        # Pipeline 2 (needed by TL pipelines)
        if (run_tl || run_emb) {
          tagList(
            hr(),
            div(class = "pipeline-block",
              tags$p(tags$strong("Step 2 — Extract global covariates from GEE"),
                tags$span(class = "badge bg-secondary ms-2", "~60 min")),
              p(class = "text-muted small",
                "Run once. Re-running is safe — completed batches are skipped."),
              tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_preanalysis.R",
  store  = "_targets_preanalysis"
)')
            )
          )
        },

        # Pipeline 3
        if (run_tl) {
          tagList(
            hr(),
            div(class = "pipeline-block",
              tags$p(tags$strong("Step 3 — Wadoux transfer learning"),
                tags$span(class = "badge bg-secondary ms-2", "~15 min")),
              tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_transfer.R",
  store  = "_targets_transfer"
)')
            )
          )
        },

        # Pipeline 4
        if (run_emb) {
          tagList(
            hr(),
            div(class = "pipeline-block",
              tags$p(tags$strong("Step 4 — Embedding transfer learning"),
                tags$span(class = "badge bg-secondary ms-2", "~30 min")),
              tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_embedding.R",
  store  = "_targets_embedding"
)')
            )
          )
        },

        hr(),
        tags$strong("If something goes wrong:"),
        tags$pre(class = "code-block",
'targets::tar_meta() |> dplyr::filter(!is.na(error)) |>
  dplyr::select(name, error)')
      )
    )
  )
}

# ── UI helpers ─────────────────────────────────────────────────────────────
pipeline_row <- function(ns, input_id, checked, enabled, label,
                          badge, badge_cls, items) {
  cb <- if (!is.null(input_id)) {
    disabled_attr <- if (!enabled) list(disabled = NA) else list()
    do.call(
      tags$input,
      c(list(
        type    = "checkbox",
        id      = if (!is.null(ns)) ns(input_id) else input_id,
        name    = if (!is.null(ns)) ns(input_id) else input_id,
        class   = "form-check-input me-2",
        checked = if (checked && enabled) NA else NULL
      ), disabled_attr)
    )
  } else {
    tags$span(class = "me-2", style = "display:inline-block; width:18px;",
      tags$input(type = "checkbox", class = "form-check-input",
        checked = NA, disabled = NA))
  }

  div(class = paste0("pipeline-option mb-3", if (!enabled) " opacity-50"),
    div(class = "d-flex align-items-start gap-2",
      div(class = "mt-1", cb),
      div(class = "flex-grow-1",
        div(class = "d-flex align-items-center gap-2 mb-1",
          tags$strong(label),
          tags$span(class = paste("badge", badge_cls), badge)
        ),
        tags$ul(class = "mb-0 small text-muted ps-3",
          lapply(items, tags$li)
        )
      )
    )
  )
}

summary_row <- function(label, value) {
  tags$tr(
    tags$th(class = "text-muted fw-normal pe-3",
      style = "width:160px; white-space:nowrap;", label),
    tags$td(value)
  )
}
