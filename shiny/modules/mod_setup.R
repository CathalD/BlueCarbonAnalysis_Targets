mod_setup_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Tell the app about your field campaign."),
      p("These details will be written to the project configuration file.")
    ),

    bslib::card(
      bslib::card_header("Project Information"),
      bslib::card_body(
        div(class = "row g-3",
          div(class = "col-md-8",
            textInput(ns("project_name"), "Project name *",
              placeholder = "e.g. Chemainus_Estuary_2026")
          ),
          div(class = "col-md-4",
            numericInput(ns("monitoring_year"), "Monitoring year *",
              value = as.integer(format(Sys.Date(), "%Y")),
              min = 2000, max = 2100, step = 1)
          )
        ),
        textInput(ns("project_location"), "Site location *",
          placeholder = "e.g. Chemainus Estuary, British Columbia, Canada")
      )
    ),

    bslib::card(
      bslib::card_header("Ecosystem Strata"),
      bslib::card_body(
        p("Enter the ecosystem codes used in your ", code("core_locations.csv"),
          " stratum column, separated by commas."),
        textInput(ns("valid_strata"), "Stratum codes *",
          value = "IM, NM, MF",
          placeholder = "e.g. IM, NM, MF"),
        div(class = "strata-legend",
          tags$table(class = "table table-sm table-borderless mb-0",
            tags$tbody(
              tags$tr(
                tags$td(tags$strong("IM")), tags$td("Intertidal Marsh — regularly flooded by tides; dense organic soils")),
              tags$tr(
                tags$td(tags$strong("NM")), tags$td("Near-Margin — upper tidal edge; transitional zone")),
              tags$tr(
                tags$td(tags$strong("MF")), tags$td("Marine Fringe — seaward edge; closest to open water")),
              tags$tr(
                tags$td(tags$strong("SG")), tags$td("Seagrass — subtidal to intertidal meadows")),
              tags$tr(
                tags$td(tags$strong("MG")), tags$td("Mangrove — woody tidal forest"))
            )
          )
        )
      )
    ),

    bslib::card(
      bslib::card_header("Google Earth Engine"),
      bslib::card_body(
        textInput(ns("gee_project"), "GEE Cloud Project ID",
          placeholder = "e.g. my-project-470316"),
        helpText(
          "Leave blank if you are only running the basic analysis (Pipeline 1). ",
          "This can be added later in ", code("blue_carbon_config.R"), "."
        )
      )
    ),

    uiOutput(ns("validation_msg"))
  )
}

mod_setup_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    output$validation_msg <- renderUI({
      v <- validate_setup()
      if (length(v$errors) > 0) {
        div(class = "alert alert-warning mt-2",
          tags$ul(lapply(v$errors, tags$li))
        )
      }
    })

    validate_setup <- reactive({
      errors <- character(0)
      if (is.null(input$project_name) || nchar(trimws(input$project_name)) == 0)
        errors <- c(errors, "Project name is required.")
      if (is.null(input$project_location) || nchar(trimws(input$project_location)) == 0)
        errors <- c(errors, "Site location is required.")
      if (length(parse_strata(input$valid_strata)) == 0)
        errors <- c(errors, "At least one stratum code is required.")
      list(errors = errors)
    })

    reactive({
      strata <- parse_strata(input$valid_strata)
      list(
        project_name     = trimws(input$project_name     %||% ""),
        project_location = trimws(input$project_location %||% ""),
        monitoring_year  = input$monitoring_year %||% as.integer(format(Sys.Date(), "%Y")),
        valid_strata     = strata,
        gee_project      = trimws(input$gee_project      %||% ""),
        ready            = length(validate_setup()$errors) == 0
      )
    })
  })
}

parse_strata <- function(s) {
  if (is.null(s) || nchar(trimws(s)) == 0) return(character(0))
  parts <- unlist(strsplit(s, "[,;[:space:]]+"))
  parts <- trimws(parts)
  parts[nchar(parts) > 0]
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
