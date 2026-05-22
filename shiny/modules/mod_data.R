mod_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Upload your field data CSV files."),
      p("One file lists your core locations; the other lists the lab measurements at each depth.")
    ),

    bslib::layout_columns(
      col_widths = c(6, 6),

      bslib::card(
        bslib::card_header("core_locations.csv"),
        bslib::card_body(
          fileInput(ns("locations_file"), NULL,
            accept   = ".csv",
            buttonLabel = "Browse…",
            placeholder = "core_locations.csv"),
          uiOutput(ns("loc_status")),
          div(style = "overflow-x: auto;",
            DT::dataTableOutput(ns("loc_preview")))
        )
      ),

      bslib::card(
        bslib::card_header("core_samples.csv"),
        bslib::card_body(
          fileInput(ns("samples_file"), NULL,
            accept   = ".csv",
            buttonLabel = "Browse…",
            placeholder = "core_samples.csv"),
          uiOutput(ns("samp_status")),
          div(style = "overflow-x: auto;",
            DT::dataTableOutput(ns("samp_preview")))
        )
      )
    ),

    bslib::card(
      bslib::card_header("CSV format reference"),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(6, 6),
          div(
            tags$strong("core_locations.csv"),
            tags$pre(class = "csv-example",
"core_id,latitude,longitude,stratum
SITE_01,48.899,-123.671,IM
SITE_02,48.901,-123.675,NM
SITE_03,48.897,-123.668,MF")
          ),
          div(
            tags$strong("core_samples.csv"),
            tags$pre(class = "csv-example",
"core_id,depth_top_cm,depth_bottom_cm,soc_g_kg,bulk_density_g_cm3
SITE_01,0,10,87.4,0.42
SITE_01,10,20,64.2,0.51
SITE_01,20,50,45.1,0.63
SITE_01,50,100,28.7,")
          )
        )
      )
    )
  )
}

mod_data_server <- function(id, setup_state) {
  moduleServer(id, function(input, output, session) {

    loc_df <- reactive({
      req(input$locations_file)
      tryCatch(
        readr::read_csv(input$locations_file$datapath, show_col_types = FALSE),
        error = function(e) NULL
      )
    })

    samp_df <- reactive({
      req(input$samples_file)
      tryCatch(
        readr::read_csv(input$samples_file$datapath, show_col_types = FALSE),
        error = function(e) NULL
      )
    })

    loc_validation <- reactive({
      df <- loc_df()
      if (is.null(df)) return(list(errors = "Could not read file.", warnings = character(0)))
      validate_locations_csv(df, valid_strata = setup_state()$valid_strata)
    })

    samp_validation <- reactive({
      df <- samp_df()
      if (is.null(df)) return(list(errors = "Could not read file.", warnings = character(0)))
      validate_samples_csv(df)
    })

    output$loc_status <- renderUI({
      df <- loc_df()
      if (is.null(df)) return(NULL)
      validation_alert(loc_validation(), nrow(df), "core(s)")
    })

    output$samp_status <- renderUI({
      df <- samp_df()
      if (is.null(df)) return(NULL)
      validation_alert(samp_validation(), nrow(df), "sample(s)")
    })

    output$loc_preview <- DT::renderDataTable({
      req(loc_df())
      DT::datatable(loc_df(),
        options = list(pageLength = 5, scrollX = TRUE, dom = "tp"),
        rownames = FALSE)
    })

    output$samp_preview <- DT::renderDataTable({
      req(samp_df())
      DT::datatable(samp_df(),
        options = list(pageLength = 5, scrollX = TRUE, dom = "tp"),
        rownames = FALSE)
    })

    reactive({
      loc  <- loc_df()
      samp <- samp_df()
      loc_ok  <- !is.null(loc)  && length(loc_validation()$errors)  == 0
      samp_ok <- !is.null(samp) && length(samp_validation()$errors) == 0
      list(
        locations = loc,
        samples   = samp,
        ready     = loc_ok && samp_ok
      )
    })
  })
}

validation_alert <- function(v, n_rows, unit) {
  has_errors   <- length(v$errors) > 0
  has_warnings <- length(v$warnings) > 0

  if (has_errors) {
    div(class = "alert alert-danger mt-2 mb-1",
      tags$strong(paste0("❌ ", length(v$errors), " error(s) found")),
      tags$ul(lapply(v$errors, tags$li))
    )
  } else if (has_warnings) {
    tagList(
      div(class = "alert alert-success mt-2 mb-1",
        paste0("✓ ", n_rows, " ", unit, " loaded successfully")
      ),
      div(class = "alert alert-warning mb-1",
        tags$strong(paste0("⚠ ", length(v$warnings), " warning(s)")),
        tags$ul(lapply(v$warnings, tags$li))
      )
    )
  } else {
    div(class = "alert alert-success mt-2 mb-1",
      paste0("✓ ", n_rows, " ", unit, " loaded — no issues found")
    )
  }
}
