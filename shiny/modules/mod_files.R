mod_files_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "d-flex justify-content-between align-items-center mb-3",
      div(
        h4(class = "mb-0", "Project files"),
        p(class = "text-muted mb-0 small",
          "All inputs and outputs for this project, and where to find them.")
      ),
      actionButton(ns("refresh"), "↺ Refresh",
        class = "btn btn-outline-secondary btn-sm")
    ),
    uiOutput(ns("content"))
  )
}

mod_files_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    refresh_count <- reactiveVal(0L)
    observeEvent(input$refresh, refresh_count(refresh_count() + 1L))

    output$content <- renderUI({
      refresh_count()

      cfg        <- read_project_config_files(project_root)
      covar_path <- cfg$COVARIATE_RASTER %||% ""

      data_raw <- file.path(project_root, "Pre-Analysis Data Preparation", "data_raw")
      locs_path <- file.path(data_raw, "core_locations.csv")
      samp_path <- file.path(data_raw, "core_samples.csv")
      aoi_path  <- file.path(data_raw, "aoi_boundary.geojson")

      tagList(

        # ── Input files ───────────────────────────────────────────────────────
        bslib::card(
          class = "mb-3",
          bslib::card_header("Input files"),
          bslib::card_body(class = "p-0",
            tags$table(class = "table table-sm table-hover mb-0",
              tags$thead(class = "table-light",
                tags$tr(
                  tags$th("File", style = "width: 45%"),
                  tags$th("Status"),
                  tags$th("Details")
                )
              ),
              tags$tbody(
                file_status_row(
                  label = "blue_carbon_config.R",
                  path  = file.path(project_root, "blue_carbon_config.R"),
                  note  = if (!is.null(cfg$PROJECT_NAME))
                    paste0("Project: ", cfg$PROJECT_NAME)
                ),
                file_status_row(
                  label = "data_raw/core_locations.csv",
                  path  = locs_path,
                  note  = csv_summary(locs_path, "location")
                ),
                file_status_row(
                  label = "data_raw/core_samples.csv",
                  path  = samp_path,
                  note  = csv_summary(samp_path, "depth measurement")
                ),
                if (nchar(trimws(covar_path)) > 0)
                  file_status_row(
                    label = paste0("covariates/", basename(covar_path)),
                    path  = covar_path,
                    note  = raster_summary(covar_path)
                  ),
                if (file.exists(aoi_path))
                  file_status_row(
                    label = "data_raw/aoi_boundary.geojson",
                    path  = aoi_path
                  )
              )
            )
          )
        ),

        # ── Reports ───────────────────────────────────────────────────────────
        bslib::card(
          class = "mb-3",
          bslib::card_header("Reports"),
          bslib::card_body(class = "p-0",
            p(class = "px-3 pt-3 pb-1 text-muted small mb-0",
              "Open HTML reports in any web browser for full results, ",
              "plots, and interpretation."),
            tags$table(class = "table table-sm table-hover mb-0",
              tags$thead(class = "table-light",
                tags$tr(
                  tags$th("Report", style = "width: 55%"),
                  tags$th("Status"),
                  tags$th("Action")
                )
              ),
              tags$tbody(
                report_row(ns, project_root, "reports/step1_nonspatial.html",
                  "Non-spatial analysis (Step 1)", "dl_s1"),
                report_row(ns, project_root, "reports/step3_random_forest.html",
                  "RF spatial maps (Pipeline 2)", "dl_s3"),
                report_row(ns, project_root, "reports/step4_transfer_learning.html",
                  "Wadoux transfer learning (Pipeline 4)", "dl_s4"),
                report_row(ns, project_root, "reports/step5_embedding_tl.html",
                  "Embedding TL (Pipeline 5)", "dl_s5")
              )
            )
          )
        ),

        # ── Spatial outputs ───────────────────────────────────────────────────
        bslib::card(
          class = "mb-3",
          bslib::card_header("Spatial outputs (GeoTIFFs)"),
          bslib::card_body(class = "p-0",
            p(class = "px-3 pt-3 pb-1 text-muted small mb-0",
              "Prediction rasters can be opened in QGIS, ArcGIS, or R with terra."),
            tags$table(class = "table table-sm table-hover mb-0",
              tags$thead(class = "table-light",
                tags$tr(tags$th("Folder", style = "width: 55%"), tags$th("Contents"))
              ),
              tags$tbody(
                raster_dir_row(project_root, "outputs/rf",
                  "RF prediction maps — one band per VM0033 depth interval + total"),
                raster_dir_row(project_root, "outputs/transfer",
                  "Wadoux TL maps — Global Prior, Transfer Final, Local Only, Difference"),
                raster_dir_row(project_root, "outputs/embedding",
                  "Embedding TL maps + 64-band AOI embedding raster")
              )
            )
          )
        ),

        # ── Folder map ────────────────────────────────────────────────────────
        bslib::card(
          bslib::card_header("Folder map — where to find everything"),
          bslib::card_body(
            tags$pre(class = "bg-light p-3 rounded mb-0", style = "font-size: 12px;",
'BlueCarbonAnalysis_Targets/
├── blue_carbon_config.R          ← your project settings
├── reports/
│   ├── step1_nonspatial.html     ← non-spatial analysis report
│   ├── step3_random_forest.html  ← RF spatial map report
│   ├── step4_transfer_learning.html
│   └── step5_embedding_tl.html
├── outputs/
│   ├── rf/                       ← RF prediction rasters (.tif)
│   ├── transfer/                 ← Wadoux TL rasters (.tif)
│   └── embedding/                ← Embedding TL rasters (.tif)
└── Pre-Analysis Data Preparation/
    ├── data_raw/
    │   ├── core_locations.csv    ← your uploaded field core locations
    │   └── core_samples.csv      ← your uploaded depth measurements
    └── covariates/
        └── YourRaster.tif        ← place covariate raster here')
          )
        )
      )
    })

    # ── Download handlers ──────────────────────────────────────────────────────
    make_dl <- function(output_id, filename) {
      output[[output_id]] <- downloadHandler(
        filename = filename,
        content  = function(file) {
          file.copy(file.path(project_root, "reports", filename), file)
        }
      )
    }
    make_dl("dl_s1", "step1_nonspatial.html")
    make_dl("dl_s3", "step3_random_forest.html")
    make_dl("dl_s4", "step4_transfer_learning.html")
    make_dl("dl_s5", "step5_embedding_tl.html")
  })
}

# ── Helpers ────────────────────────────────────────────────────────────────────

read_project_config_files <- function(project_root) {
  config_path <- file.path(project_root, "blue_carbon_config.R")
  if (!file.exists(config_path)) return(list())
  e <- new.env(parent = emptyenv())
  tryCatch(
    { source(config_path, local = e, echo = FALSE); as.list(e) },
    error = function(x) list()
  )
}

file_status_row <- function(label, path, note = NULL) {
  exists <- file.exists(path)
  fi     <- if (exists) file.info(path) else NULL
  detail <- if (!is.null(note) && nchar(trimws(note)) > 0) {
    note
  } else if (exists) {
    paste0(
      round(fi$size / 1024, 1), " KB",
      " — ", format(fi$mtime, "%Y-%m-%d")
    )
  } else {
    "Not found"
  }
  tags$tr(
    tags$td(tags$code(label, style = "font-size:12px;")),
    tags$td(
      if (exists) span(class = "text-success fw-bold", "✓")
      else span(class = "text-muted", "—")
    ),
    tags$td(class = "text-muted small", detail)
  )
}

csv_summary <- function(path, unit) {
  if (!file.exists(path)) return(NULL)
  n <- tryCatch(
    nrow(readr::read_csv(path, show_col_types = FALSE)),
    error = function(e) NULL
  )
  if (is.null(n)) return(NULL)
  paste0(n, " ", unit, if (n != 1) "s")
}

raster_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    r <- terra::rast(path)
    paste0(terra::nlyr(r), " bands, ~",
           round(mean(terra::res(r)), 0), " m resolution")
  }, error = function(e) NULL)
}

report_row <- function(ns, project_root, rel_path, label, dl_id) {
  full_path <- file.path(project_root, rel_path)
  exists    <- file.exists(full_path)
  tags$tr(
    tags$td(label),
    tags$td(
      if (exists) span(class = "text-success fw-bold", "✓ Ready")
      else span(class = "text-muted small", "Not yet created")
    ),
    tags$td(
      if (exists)
        downloadButton(ns(dl_id), "Download",
          class = "btn btn-outline-primary btn-sm py-0")
      else
        span(class = "text-muted small", "Run the pipeline first")
    )
  )
}

raster_dir_row <- function(project_root, rel_path, description) {
  full_path <- file.path(project_root, rel_path)
  n_tifs    <- if (dir.exists(full_path)) {
    length(list.files(full_path, pattern = "\\.tif$", recursive = TRUE))
  } else 0L

  contents <- if (n_tifs > 0) {
    paste0(n_tifs, " .tif file", if (n_tifs > 1) "s")
  } else {
    span(class = "text-muted small", "Empty — run the pipeline to populate")
  }

  tags$tr(
    tags$td(tags$code(rel_path, style = "font-size: 12px;")),
    tags$td(class = "small",
      div(contents),
      div(class = "text-muted", description)
    )
  )
}
