suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve project root — shiny/ sits one level below it
PROJECT_ROOT <- normalizePath(file.path(getwd(), ".."))

# Source helpers and modules
source(file.path(getwd(), "R",       "validate_csv.R"))
source(file.path(getwd(), "R",       "write_config.R"))
source(file.path(getwd(), "modules", "mod_setup.R"))
source(file.path(getwd(), "modules", "mod_data.R"))
source(file.path(getwd(), "modules", "mod_raster.R"))
source(file.path(getwd(), "modules", "mod_finish.R"))
source(file.path(getwd(), "modules", "mod_guide.R"))
source(file.path(getwd(), "modules", "mod_files.R"))

# ── Step definitions (Setup wizard) ─────────────────────────────────────────
STEPS <- list(
  list(id = "step1", label = "1. Project Setup"),
  list(id = "step2", label = "2. Field Data"),
  list(id = "step3", label = "3. Raster & AOI"),
  list(id = "step4", label = "4. Save & Configure")
)

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = div(
    class = "d-flex align-items-center gap-2",
    tags$svg(
      viewBox = "0 0 32 32", width = "26", height = "26",
      tags$circle(cx = "16", cy = "16", r = "16", fill = "white"),
      tags$text(
        x = "16", y = "22", `text-anchor` = "middle",
        fill = "#2c7a4b", `font-size` = "14", `font-weight` = "bold", "BC"
      )
    ),
    span("Blue Carbon Analysis", style = "color: white; font-weight: 600;")
  ),
  theme    = bs_theme(bootswatch = "flatly", version = 5),
  bg       = "#2c7a4b",
  inverse  = TRUE,
  includeCSS(file.path(getwd(), "www", "custom.css")),

  # ── Tab 1: Setup wizard ──────────────────────────────────────────────────
  nav_panel(
    title = "Setup",
    icon  = icon("sliders"),

    uiOutput("step_indicator"),

    div(class = "container-fluid px-4 py-4",
      tabsetPanel(
        id   = "wizard_tabs",
        type = "hidden",
        tabPanelBody("step1", mod_setup_ui("setup")),
        tabPanelBody("step2", mod_data_ui("data")),
        tabPanelBody("step3", mod_raster_ui("raster")),
        tabPanelBody("step4", mod_finish_ui("finish"))
      ),
      div(class = "wizard-nav d-flex justify-content-between mt-4",
        uiOutput("btn_back_ui"),
        uiOutput("btn_next_ui")
      )
    )
  ),

  # ── Tab 2: Run guide ─────────────────────────────────────────────────────
  nav_panel(
    title = "Run",
    icon  = icon("terminal"),
    div(class = "container-fluid px-4 py-4",
      mod_guide_ui("guide")
    )
  ),

  # ── Tab 3: Project files ──────────────────────────────────────────────────
  nav_panel(
    title = "Outputs",
    icon  = icon("folder-open"),
    div(class = "container-fluid px-4 py-4",
      mod_files_ui("files")
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  current_step <- reactiveVal(1L)

  # Module servers
  setup_state  <- mod_setup_server("setup")
  data_state   <- mod_data_server("data", setup_state)
  raster_state <- mod_raster_server("raster", PROJECT_ROOT)
  mod_finish_server("finish", setup_state, data_state, raster_state, PROJECT_ROOT)
  mod_guide_server("guide", PROJECT_ROOT)
  mod_files_server("files", PROJECT_ROOT)

  # ── Step indicator (inside Setup tab) ────────────────────────────────────
  output$step_indicator <- renderUI({
    step  <- current_step()
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

  # Sync tab panel with step counter
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

  # Next / Back
  observeEvent(input$btn_next, {
    step <- current_step()
    if (step < length(STEPS)) current_step(step + 1L)
  })

  observeEvent(input$btn_back, {
    step <- current_step()
    if (step > 1L) current_step(step - 1L)
  })

  output$btn_back_ui <- renderUI({
    if (current_step() <= 1L) return(div())
    actionButton("btn_back", "← Back", class = "btn btn-outline-secondary")
  })

  output$btn_next_ui <- renderUI({
    step  <- current_step()
    if (step >= length(STEPS)) return(div())
    ready <- step_ready()
    cls   <- if (ready) "btn btn-primary" else "btn btn-primary disabled"
    tip   <- if (!ready) switch(step,
      `1` = "Complete project name and location to continue.",
      `2` = "Upload both CSV files without errors to continue.",
      `3` = "Confirm raster setting to continue.",
      NULL
    )
    tagList(
      if (!is.null(tip)) p(class = "text-muted small mb-1 text-end", tip),
      actionButton("btn_next", "Next →", class = cls)
    )
  })
}

shinyApp(ui, server)
