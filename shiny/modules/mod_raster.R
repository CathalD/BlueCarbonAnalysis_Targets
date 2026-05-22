mod_raster_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Satellite data — optional."),
      p("The basic non-spatial analysis (depth harmonization, stratum summaries, ",
        "VM0033 carbon stocks) ", tags$strong("does not need a raster."),
        " Come back to this step later if you want spatial prediction maps.")
    ),

    bslib::card(
      bslib::card_header("Covariate Raster — for spatial prediction maps"),
      bslib::card_body(
        checkboxInput(ns("skip_raster"),
          "Skip — run non-spatial analysis only (recommended if you don't have the raster yet)",
          value = TRUE),

        conditionalPanel(
          condition = paste0("!input['", ns("skip_raster"), "']"),
          div(class = "mt-3",
            p("Place your GEE covariate .tif in ",
              code("Pre-Analysis Data Preparation/covariates/"),
              " then enter the filename below."),
            div(class = "row g-2 align-items-end",
              div(class = "col",
                textInput(ns("raster_name"), NULL,
                  placeholder = "e.g. BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif")
              ),
              div(class = "col-auto",
                actionButton(ns("check_raster"), "Check file",
                  class = "btn btn-secondary btn-sm")
              )
            ),
            uiOutput(ns("raster_status"))
          )
        )
      )
    ),

    bslib::card(
      bslib::card_header("Site Boundary — for area-weighted totals"),
      bslib::card_body(
        checkboxInput(ns("skip_aoi"),
          "Skip — return per-stratum carbon densities only (kg C/m²)",
          value = TRUE),
        conditionalPanel(
          condition = paste0("!input['", ns("skip_aoi"), "']"),
          div(class = "mt-3",
            fileInput(ns("aoi_file"), NULL,
              accept      = c(".geojson", ".json", ".gpkg", ".zip"),
              buttonLabel = "Upload boundary file",
              placeholder = "GeoJSON / GPKG / zipped shapefile")
          )
        ),
        uiOutput(ns("aoi_status"))
      )
    )
  )
}

mod_raster_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {

    raster_info <- reactiveVal(NULL)

    observeEvent(input$check_raster, {
      rname <- trimws(input$raster_name %||% "")
      if (nchar(rname) == 0) {
        raster_info(list(found = FALSE,
          msg = "Enter a filename or path above, then click Check file."))
        return()
      }

      candidates <- c(
        rname,
        file.path(project_root, rname),
        file.path(project_root, "Pre-Analysis Data Preparation", "covariates", rname)
      )
      found_path <- Find(file.exists, candidates)

      if (is.null(found_path)) {
        raster_info(list(found = FALSE, msg = paste0(
          "File not found. Make sure the .tif is in:\n",
          "  Pre-Analysis Data Preparation/covariates/\n",
          "and the filename matches exactly."
        )))
        return()
      }

      info <- tryCatch({
        suppressPackageStartupMessages(library(terra))
        r <- terra::rast(found_path)
        list(
          found  = TRUE,
          path   = found_path,
          bands  = terra::nlyr(r),
          crs    = tryCatch(terra::crs(r, describe = TRUE)$name, error = function(e) "unknown"),
          res_m  = round(mean(terra::res(r)), 1),
          msg    = paste0("Found: ", basename(found_path))
        )
      }, error = function(e) {
        list(found = FALSE, msg = paste("Could not read raster:", conditionMessage(e)))
      })

      raster_info(info)
    })

    output$raster_status <- renderUI({
      info <- raster_info()
      if (is.null(info)) return(helpText("Click 'Check file' to validate."))
      if (!info$found) {
        div(class = "alert alert-danger mt-2",
          tags$strong("❌ File not found"), tags$br(),
          tags$small(tags$pre(info$msg))
        )
      } else {
        div(class = "alert alert-success mt-2",
          tags$strong(paste0("✓ ", info$msg)), tags$br(),
          paste0(info$bands, " bands | CRS: ", info$crs,
                 " | Resolution: ~", info$res_m, " m")
        )
      }
    })

    output$aoi_status <- renderUI({
      if (isTRUE(input$skip_aoi)) return(NULL)
      if (!is.null(input$aoi_file))
        div(class = "alert alert-success mt-2",
          paste0("✓ ", input$aoi_file$name, " uploaded"))
    })

    reactive({
      skip_r <- isTRUE(input$skip_raster)
      info    <- raster_info()

      raster_path <- if (!skip_r && !is.null(info) && isTRUE(info$found)) {
        info$path
      } else ""

      aoi_source <- if (!isTRUE(input$skip_aoi) && !is.null(input$aoi_file)) {
        input$aoi_file$datapath
      } else NULL

      aoi_dest <- if (!is.null(aoi_source)) {
        file.path(project_root, "Pre-Analysis Data Preparation", "data_raw",
                  "aoi_boundary.geojson")
      } else NULL

      list(
        raster_path  = raster_path,
        has_raster   = nchar(raster_path) > 0,
        aoi_source   = aoi_source,
        aoi_dest     = aoi_dest,
        ready        = TRUE   # raster is always optional
      )
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
