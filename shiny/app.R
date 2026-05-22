suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(readr)
})

# Resolve project root — shiny/ sits one level below it
PROJECT_ROOT <- normalizePath(file.path(getwd(), ".."))

# Source helpers and modules
source(file.path(getwd(), "R",       "validate_csv.R"))
source(file.path(getwd(), "R",       "write_config.R"))
source(file.path(getwd(), "modules", "mod_setup.R"))
source(file.path(getwd(), "modules", "mod_data.R"))
source(file.path(getwd(), "modules", "mod_raster.R"))
source(file.path(getwd(), "modules", "mod_finish.R"))

# ── Step definitions ────────────────────────────────────────────────────────
STEPS <- list(
  list(id = "step1", label = "1. Project Setup"),
  list(id = "step2", label = "2. Field Data"),
  list(id = "step3", label = "3. Raster & AOI"),
  list(id = "step4", label = "4. Save & Run")
)

# ── UI ──────────────────────────────────────────────────────────────────────
ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly", version = 5),
  includeCSS(file.path(getwd(), "www", "custom.css")),

  # Header
  div(class = "app-header py-3 px-4",
    div(class = "d-flex align-items-center",
      div(class = "me-3",
        tags$svg(viewBox = "0 0 32 32", width = "40", height = "40",
          tags$circle(cx = "16", cy = "16", r = "16", fill = "#2c7a4b"),
          tags$text(x = "16", y = "22", `text-anchor` = "middle",
            fill = "white", `font-size` = "16", `font-weight` = "bold", "BC")
        )
      ),
      div(
        tags$h1(class = "h4 mb-0 text-white", "Blue Carbon Analysis"),
        tags$p(class = "mb-0 small", style = "color: rgba(255,255,255,0.7);",
          "Project Setup Wizard")
      )
    )
  ),

  # Step indicator
  uiOutput("step_indicator"),

  # Wizard content
  div(class = "container-fluid px-4 py-4",
    tabsetPanel(
      id   = "wizard_tabs",
      type = "hidden",
      tabPanelBody("step1", mod_setup_ui("setup")),
      tabPanelBody("step2", mod_data_ui("data")),
      tabPanelBody("step3", mod_raster_ui("raster")),
      tabPanelBody("step4", mod_finish_ui("finish"))
    ),

    # Navigation
    div(class = "wizard-nav d-flex justify-content-between mt-4",
      uiOutput("btn_back_ui"),
      uiOutput("btn_next_ui")
    )
  )
)

# ── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  current_step <- reactiveVal(1L)

  # Module servers
  setup_state  <- mod_setup_server("setup")
  data_state   <- mod_data_server("data", setup_state)
  raster_state <- mod_raster_server("raster", PROJECT_ROOT)
  mod_finish_server("finish", setup_state, data_state, raster_state, PROJECT_ROOT)

  # Step indicator
  output$step_indicator <- renderUI({
    step <- current_step()
    items <- lapply(seq_along(STEPS), function(i) {
      is_active   <- i == step
      is_complete <- i < step
      cls <- paste(c(
        "step-item",
        if (is_active)   "active",
        if (is_complete) "complete"
      ), collapse = " ")

      div(class = cls,
        div(class = "step-number", if (is_complete) "✓" else as.character(i)),
        div(class = "step-label",  STEPS[[i]]$label)
      )
    })

    div(class = "step-indicator-bar px-4 py-2",
      div(class = "d-flex gap-3 align-items-center", items)
    )
  })

  # Sync tab with step counter
  observeEvent(current_step(), {
    updateTabsetPanel(session, "wizard_tabs",
      selected = paste0("step", current_step()))
  })

  # Can we advance from the current step?
  step_ready <- reactive({
    switch(current_step(),
      `1` = setup_state()$ready,
      `2` = data_state()$ready,
      `3` = raster_state()$ready,
      `4` = TRUE,
      FALSE
    )
  })

  # Next
  observeEvent(input$btn_next, {
    step <- current_step()
    if (step < length(STEPS)) current_step(step + 1L)
  })

  # Back
  observeEvent(input$btn_back, {
    step <- current_step()
    if (step > 1L) current_step(step - 1L)
  })

  # Back button UI
  output$btn_back_ui <- renderUI({
    if (current_step() <= 1L) return(div())
    actionButton("btn_back", "← Back", class = "btn btn-outline-secondary")
  })

  # Next button UI — disabled (not hidden) when step isn't ready
  output$btn_next_ui <- renderUI({
    step <- current_step()
    if (step >= length(STEPS)) return(div())

    ready <- step_ready()
    label <- "Next →"
    cls   <- if (ready) "btn btn-primary" else "btn btn-primary disabled"

    tip <- if (!ready) {
      switch(step,
        `1` = "Complete project name and location to continue.",
        `2` = "Upload both CSV files without errors to continue.",
        `3` = "Confirm raster setting to continue.",
        NULL
      )
    }

    tagList(
      if (!is.null(tip)) p(class = "text-muted small mb-1 text-end", tip),
      actionButton("btn_next", label, class = cls)
    )
  })
}

shinyApp(ui, server)
