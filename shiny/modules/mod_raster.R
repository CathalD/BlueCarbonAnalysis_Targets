mod_raster_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Point the app to your satellite covariate raster."),
      p("The covariate raster is a multi-band .tif file exported from Google Earth Engine.",
        "Place it in ", code("Pre-Analysis Data Preparation/covariates/"),
        " before continuing.")
    ),

    bslib::card(
      bslib::card_header("Covariate Raster"),
      bslib::card_body(
        div(class = "row g-2 align-items-end",
          div(class = "col",
            textInput(ns("raster_name"), "Raster filename or full path",
              placeholder = "e.g. BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif")
          ),
          div(class = "col-auto",
            actionButton(ns("check_raster"), "Check file",
              class = "btn btn-secondary", icon = icon("search"))
          )
        ),
        helpText("If the file is in Pre-Analysis Data Preparation/covariates/, ",
                 "just enter the filename. Otherwise paste the full path."),
        uiOutput(ns("raster_status")),
        checkboxInput(ns("skip_raster"), "Skip — I don't have a raster yet (basic stats only)", value = FALSE)
      )
    ),

    bslib::card(
      bslib::card_header("Site Boundary (optional)"),
      bslib::card_body(
        p("Upload a site boundary file to get total carbon stock estimates for the whole site.",
          "Without it, the pipeline returns per-stratum carbon densities only."),
        fileInput(ns("aoi_file"), NULL,
          accept      = c(".geojson", ".json", ".gpkg", ".zip"),
          buttonLabel = "Upload boundary",
          placeholder = "GeoJSON / GPKG / zipped shapefile"),
        checkboxInput(ns("skip_aoi"), "Skip — I don't have a boundary file", value = TRUE),
        uiOutput(ns("aoi_status"))
      )
    )
  )
}

mod_raster_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {

    raster_info <- reactiveVal(NULL)

    observeEvent(input$check_raster, {
      rname <- trimws(input$raster_name)
      if (nchar(rname) == 0) {
        raster_info(list(found = FALSE, msg = "Enter a filename or path above, then click Check file."))
        return()
      }

      candidates <- c(
        rname,
        file.path(project_root, rname),
        file.path(project_root, "Pre-Analysis Data Preparation", "covariates", rname)
      )
      found_path <- Find(file.exists, candidates)

      if (is.null(found_path)) {
        raster_info(list(
          found = FALSE,
          msg   = paste0(
            "File not found. Make sure the .tif is in:\n",
            "  Pre-Analysis Data Preparation/covariates/\n",
            "and the filename matches exactly (including capitalisation)."
          )
        ))
        return()
      }

      info <- tryCatch({
        suppressPackageStartupMessages(library(terra))
        r <- terra::rast(found_path)
        list(
          found   = TRUE,
          path    = found_path,
          bands   = terra::nlyr(r),
          crs     = tryCatch(terra::crs(r, describe = TRUE)$name, error = function(e) "unknown"),
          res_m   = round(mean(terra::res(r)), 1),
          msg     = paste0("Found: ", basename(found_path))
        )
      }, error = function(e) {
        list(found = FALSE, msg = paste("File found but could not be read:", conditionMessage(e)))
      })

      raster_info(info)
    })

    output$raster_status <- renderUI({
      if (isTRUE(input$skip_raster)) {
        return(div(class = "alert alert-warning mt-2",
          "⚠ Skipping raster — pipeline will produce depth summaries and plots only. ",
          "Spatial prediction maps (Step 3 RF, Steps 4–5 TL) will not be generated."
        ))
      }

      info <- raster_info()
      if (is.null(info)) return(helpText("Click 'Check file' to validate."))

      if (!info$found) {
        div(class = "alert alert-danger mt-2",
          tags$strong("❌ File not found"),
          tags$br(),
          tags$small(pre(info$msg))
        )
      } else {
        div(class = "alert alert-success mt-2",
          tags$strong(paste0("✓ ", info$msg)),
          tags$br(),
          paste0(info$bands, " bands | CRS: ", info$crs, " | Resolution: ~", info$res_m, " m")
        )
      }
    })

    output$aoi_status <- renderUI({
      if (isTRUE(input$skip_aoi)) {
        div(class = "alert alert-info mt-2",
          "ℹ No boundary file — the pipeline will return per-stratum carbon ",
          "densities (kg C/m²) without an area-weighted site total."
        )
      } else if (!is.null(input$aoi_file)) {
        div(class = "alert alert-success mt-2",
          paste0("✓ Boundary uploaded: ", input$aoi_file$name)
        )
      }
    })

    reactive({
      info <- raster_info()

      raster_path <- if (isTRUE(input$skip_raster)) {
        ""
      } else if (!is.null(info) && isTRUE(info$found)) {
        info$path
      } else {
        ""
      }

      aoi_source <- if (!isTRUE(input$skip_aoi) && !is.null(input$aoi_file)) {
        input$aoi_file$datapath
      } else NULL

      aoi_dest <- if (!is.null(aoi_source)) {
        file.path(project_root, "Pre-Analysis Data Preparation", "data_raw", "aoi_boundary.geojson")
      } else NULL

      raster_ok <- isTRUE(input$skip_raster) || (nchar(raster_path) > 0)

      list(
        raster_path = raster_path,
        aoi_source  = aoi_source,
        aoi_dest    = aoi_dest,
        ready       = raster_ok
      )
    })
  })
}
