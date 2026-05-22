suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(readr)
  library(callr)
})

# Resolve project root: shiny/ is one level below the project root
PROJECT_ROOT <- normalizePath(file.path(getwd(), ".."))

# Source helpers and modules
source(file.path(getwd(), "R",       "validate_csv.R"))
source(file.path(getwd(), "R",       "write_config.R"))
source(file.path(getwd(), "modules", "mod_setup.R"))
source(file.path(getwd(), "modules", "mod_data.R"))
source(file.path(getwd(), "modules", "mod_raster.R"))
source(file.path(getwd(), "modules", "mod_run.R"))
source(file.path(getwd(), "modules", "mod_results.R"))

# ── Step metadata ──────────────────────────────────────────────────────────────
STEPS <- list(
  list(id = "step1", label = "1. Project Setup",  icon = "gear"),
  list(id = "step2", label = "2. Field Data",      icon = "upload"),
  list(id = "step3", label = "3. Raster & AOI",    icon = "map"),
  list(id = "step4", label = "4. Run Pipeline 1",  icon = "play-circle"),
  list(id = "step5", label = "5. Results",          icon = "bar-chart")
)

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly", version = 5),
  includeCSS(file.path(getwd(), "www", "custom.css")),

  # Header
  div(class = "app-header py-3 px-4 mb-0",
    div(class = "d-flex align-items-center",
      div(class = "me-3",
        tags$svg(class = "app-logo", viewBox = "0 0 32 32", width = "40", height = "40",
          tags$circle(cx = "16", cy = "16", r = "16", fill = "#2c7a4b"),
          tags$text(x = "16", y = "22", `text-anchor` = "middle",
            fill = "white", `font-size` = "18", `font-weight` = "bold", "BC")
        )
      ),
      div(
        tags$h1(class = "h4 mb-0 text-white", "Blue Carbon Analysis"),
        tags$p(class = "mb-0 text-white-50 small", "Project Setup Wizard")
      )
    )
  ),

  # Step indicator
  uiOutput("step_indicator"),

  # Main content
  div(class = "container-fluid px-4 py-3",
    tabsetPanel(
      id   = "wizard_tabs",
      type = "hidden",
      tabPanelBody("step1", mod_setup_ui("setup")),
      tabPanelBody("step2", mod_data_ui("data")),
      tabPanelBody("step3", mod_raster_ui("raster")),
      tabPanelBody("step4", mod_run_ui("run")),
      tabPanelBody("step5", mod_results_ui("results"))
    ),

    # Navigation buttons
    div(class = "wizard-nav d-flex justify-content-between mt-4",
      uiOutput("btn_back_ui"),
      uiOutput("btn_next_ui")
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  current_step <- reactiveVal(1L)
  save_error   <- reactiveVal(NULL)

  # Module servers
  setup_state  <- mod_setup_server("setup")
  data_state   <- mod_data_server("data", setup_state)
  raster_state <- mod_raster_server("raster", PROJECT_ROOT)
  run_state    <- mod_run_server("run", PROJECT_ROOT)
  mod_results_server("results", PROJECT_ROOT)

  # Step indicator
  output$step_indicator <- renderUI({
    step <- current_step()
    items <- lapply(seq_along(STEPS), function(i) {
      s <- STEPS[[i]]
      is_active   <- i == step
      is_complete <- i < step
      cls <- paste0(
        "step-item",
        if (is_active)   " active",
        if (is_complete) " complete"
      )
      div(class = cls,
        div(class = "step-number",
          if (is_complete) "✓" else as.character(i)
        ),
        div(class = "step-label", s$label)
      )
    })

    div(class = "step-indicator-bar px-4 py-2",
      div(class = "d-flex justify-content-between align-items-center",
        items
      )
    )
  })

  # Sync tab panel with current step
  observeEvent(current_step(), {
    updateTabsetPanel(session, "wizard_tabs",
      selected = paste0("step", current_step()))
  })

  # ── Validation per step ──
  step_ready <- reactive({
    switch(current_step(),
      `1` = setup_state()$ready,
      `2` = data_state()$ready,
      `3` = raster_state()$ready,
      `4` = run_state()$done,
      `5` = TRUE,
      FALSE
    )
  })

  # ── Save project files (triggered on step 3 → 4) ──
  save_project <- function() {
    s <- setup_state()
    d <- data_state()
    r <- raster_state()

    tryCatch({
      data_raw <- file.path(PROJECT_ROOT, "Pre-Analysis Data Preparation", "data_raw")
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
        output_path      = file.path(PROJECT_ROOT, "blue_carbon_config.R")
      )

      save_error(NULL)
      TRUE
    }, error = function(e) {
      save_error(conditionMessage(e))
      FALSE
    })
  }

  # ── Next button ──
  observeEvent(input$btn_next, {
    step <- current_step()

    if (step == 3L) {
      ok <- save_project()
      if (!ok) return()
    }

    if (step < 5L) current_step(step + 1L)
  })

  # ── Back button ──
  observeEvent(input$btn_back, {
    step <- current_step()
    if (step > 1L) current_step(step - 1L)
  })

  # ── Button UI ──
  output$btn_back_ui <- renderUI({
    if (current_step() <= 1L) return(div())
    actionButton("btn_back", "← Back", class = "btn btn-outline-secondary")
  })

  output$btn_next_ui <- renderUI({
    step <- current_step()
    if (step == 4L) return(div())  # advance automatically after pipeline done

    label <- if (step == 5L) "Done" else "Next →"
    cls   <- if (step_ready()) "btn btn-primary" else "btn btn-primary disabled"

    tagList(
      if (!is.null(save_error())) {
        div(class = "text-danger me-3 align-self-center small", save_error())
      },
      actionButton("btn_next", label, class = cls)
    )
  })

  # Auto-advance from step 4 → 5 when pipeline completes
  observeEvent(run_state()$done, {
    if (isTRUE(run_state()$done) && current_step() == 4L) {
      current_step(5L)
    }
  })
}

shinyApp(ui, server)
