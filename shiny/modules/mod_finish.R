mod_finish_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Review your setup and save."),
      p("Check the summary below, then click ", tags$strong("Save Setup Files"),
        " to write your configuration. After that you'll get the exact commands ",
        "to paste into RStudio to run the analysis.")
    ),

    uiOutput(ns("summary_card")),

    uiOutput(ns("save_result"))
  )
}

mod_finish_server <- function(id, setup_state, data_state, raster_state, project_root) {
  moduleServer(id, function(input, output, session) {

    saved <- reactiveVal(FALSE)
    save_error_msg <- reactiveVal(NULL)

    output$summary_card <- renderUI({
      s <- setup_state()
      d <- data_state()
      r <- raster_state()

      n_cores   <- if (!is.null(d$locations)) nrow(d$locations) else 0
      n_samples <- if (!is.null(d$samples))   nrow(d$samples)   else 0
      has_raster <- nchar(r$raster_path) > 0
      has_aoi    <- !is.null(r$aoi_source)

      bslib::card(
        bslib::card_header("Setup Summary"),
        bslib::card_body(
          tags$table(class = "table table-sm mb-0",
            tags$tbody(
              summary_row("Project name",    s$project_name),
              summary_row("Site location",   s$project_location),
              summary_row("Monitoring year", as.character(s$monitoring_year)),
              summary_row("Strata",          paste(s$valid_strata, collapse = ", ")),
              summary_row("GEE project",
                if (nchar(s$gee_project) > 0) s$gee_project
                else tags$em("Not set — required for Pipelines 2, 3, 4")),
              summary_row("Core locations",  paste0(n_cores, " cores")),
              summary_row("Sample records",  paste0(n_samples, " depth intervals")),
              summary_row("Covariate raster",
                if (has_raster) tags$span(class = "text-success", "✓ Configured")
                else tags$span(class = "text-warning",
                  "⚠ Not set — spatial maps will not be generated")),
              summary_row("Site boundary",
                if (has_aoi) tags$span(class = "text-success", "✓ Uploaded")
                else tags$span(class = "text-muted", "Not provided (optional)"))
            )
          ),
          hr(),
          tags$strong("Files that will be written:"),
          tags$ul(class = "mb-0 mt-1",
            tags$li(code("blue_carbon_config.R"), " — project configuration"),
            tags$li(code("Pre-Analysis Data Preparation/data_raw/core_locations.csv")),
            tags$li(code("Pre-Analysis Data Preparation/data_raw/core_samples.csv")),
            if (has_aoi) tags$li(code("Pre-Analysis Data Preparation/data_raw/aoi_boundary.geojson"))
          ),
          div(class = "mt-3",
            actionButton(ns("save_btn"), "Save Setup Files",
              class = "btn btn-success btn-lg",
              icon  = icon("save"))
          )
        )
      )
    })

    observeEvent(input$save_btn, {
      s <- setup_state()
      d <- data_state()
      r <- raster_state()

      result <- tryCatch({
        data_raw <- file.path(project_root,
          "Pre-Analysis Data Preparation", "data_raw")
        dir.create(data_raw, showWarnings = FALSE, recursive = TRUE)

        readr::write_csv(d$locations,
          file.path(data_raw, "core_locations.csv"))
        readr::write_csv(d$samples,
          file.path(data_raw, "core_samples.csv"))

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
      }, error = function(e) {
        save_error_msg(conditionMessage(e))
        FALSE
      })

      if (result) {
        saved(TRUE)
        save_error_msg(NULL)
      }
    })

    output$save_result <- renderUI({
      if (!is.null(save_error_msg())) {
        return(div(class = "alert alert-danger mt-3",
          tags$strong("Save failed: "), save_error_msg()
        ))
      }

      req(saved())
      s <- setup_state()

      tagList(
        div(class = "alert alert-success mt-3",
          tags$strong("✓ Setup files saved."),
          " Close this app, open the project in RStudio, and run the commands below."
        ),

        bslib::card(
          class = "mt-3 border-primary",
          bslib::card_header(class = "bg-primary text-white",
            "Run the analysis — paste these commands into the RStudio Console"),
          bslib::card_body(
            p(class = "text-muted mb-3",
              "Run these in order. Each command can be re-run safely — ",
              "completed steps are always skipped."),

            tags$strong("Step 1 — Core data processing & local RF map (5–15 min):"),
            tags$pre(class = "code-block",
"targets::tar_make()"),

            tags$strong("Check what will run before you start:"),
            tags$pre(class = "code-block",
"targets::tar_visnetwork()   # opens an interactive dependency graph"),

            if (nchar(s$gee_project) > 0) {
              tagList(
                hr(),
                div(class = "alert alert-info py-2",
                  tags$strong("Steps 2–4 require Google Earth Engine."),
                  " Run this once to authenticate (a browser window will open):"
                ),
                tags$pre(class = "code-block",
paste0('library(rgee)
ee_Initialize(user = "your.email@gmail.com", drive = TRUE)')),

                tags$strong("Step 2 — Extract global covariates from GEE (~60 min):"),
                tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_preanalysis.R",
  store  = "_targets_preanalysis"
)'),

                tags$strong("Step 3 — Transfer learning with global wetland data (~15 min):"),
                tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_transfer.R",
  store  = "_targets_transfer"
)'),

                tags$strong("Step 4 — Embedding-based transfer learning (~30 min, optional):"),
                tags$pre(class = "code-block",
'targets::tar_make(
  script = "_targets_embedding.R",
  store  = "_targets_embedding"
)')
              )
            } else {
              tagList(
                hr(),
                div(class = "alert alert-info py-2",
                  "To run transfer learning (Steps 2–4), add your GEE Cloud Project ID ",
                  "to ", code("blue_carbon_config.R"), " under ", code("GEE_PROJECT"),
                  " and re-run this wizard, or edit the file directly."
                )
              )
            },

            hr(),
            tags$strong("If a step fails, check what went wrong:"),
            tags$pre(class = "code-block",
'targets::tar_meta() |> dplyr::filter(!is.na(error)) |>
  dplyr::select(name, error)')
          )
        )
      )
    })
  })
}

summary_row <- function(label, value) {
  tags$tr(
    tags$th(class = "text-muted fw-normal", style = "width: 200px;", label),
    tags$td(value)
  )
}
