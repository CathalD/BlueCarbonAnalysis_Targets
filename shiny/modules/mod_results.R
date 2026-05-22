mod_results_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Pipeline 1 is complete."),
      p("Your project is now configured and the core data processing has finished. ",
        "Review the results below, then move to the advanced pipelines using the ",
        "copy-paste commands at the bottom of this page.")
    ),

    bslib::card(
      bslib::card_body(class = "p-0",
        tabsetPanel(
          tabPanel("Stratum Summary",
            div(class = "p-3",
              p(class = "text-muted",
                "Mean carbon stock (kg C/m²) at each VM0033 depth interval, ",
                "grouped by ecosystem stratum. These are the core VM0033 reporting values."),
              DT::dataTableOutput(ns("summary_table"))
            )
          ),
          tabPanel("Depth Profiles",
            div(class = "p-3",
              p(class = "text-muted",
                "SOC vs. depth profiles for each core. Carbon should generally ",
                "decrease with depth — flag any cores where it increases strongly."),
              plotOutput(ns("depth_plot"), height = "480px")
            )
          ),
          tabPanel("Download Reports",
            div(class = "p-3",
              p("Open the full HTML reports in your browser for complete results."),
              uiOutput(ns("report_links"))
            )
          )
        )
      )
    ),

    bslib::card(
      class = "mt-4 border-primary",
      bslib::card_header(
        class = "bg-primary text-white",
        "Ready for advanced analysis — copy these commands into RStudio"
      ),
      bslib::card_body(
        p("Close this app, open the project in RStudio, and paste these commands ",
          "into the R Console to run the remaining pipelines."),

        tags$strong("One-time GEE authentication (a browser window will open):"),
        tags$pre(class = "bg-light p-3 rounded",
"library(rgee)
ee_Initialize(user = \"your.email@gmail.com\", drive = TRUE)"),

        tags$strong("Pipeline 2 — Extract global covariates from GEE (~60 min):"),
        tags$pre(class = "bg-light p-3 rounded",
"targets::tar_make(
  script = \"_targets_preanalysis.R\",
  store  = \"_targets_preanalysis\"
)"),

        tags$strong("Pipeline 3 — Wadoux transfer learning (~15 min):"),
        tags$pre(class = "bg-light p-3 rounded",
"targets::tar_make(
  script = \"_targets_transfer.R\",
  store  = \"_targets_transfer\"
)"),

        tags$strong("Pipeline 4 — Embedding transfer learning (~30 min, optional):"),
        tags$pre(class = "bg-light p-3 rounded",
"targets::tar_make(
  script = \"_targets_embedding.R\",
  store  = \"_targets_embedding\"
)"),

        div(class = "alert alert-info mt-3",
          tags$strong("ℹ Your project folder is ready."),
          " All configuration files and data have been saved. ",
          "Open ", code("BlueCarbonAnalysis_Targets.Rproj"), " in RStudio to continue."
        )
      )
    )
  )
}

mod_results_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {

    read_target <- function(name) {
      tryCatch(
        targets::tar_read(name, store = file.path(project_root, "_targets")),
        error = function(e) NULL
      )
    }

    output$summary_table <- DT::renderDataTable({
      df <- read_target("stratum_summary")
      req(!is.null(df))

      DT::datatable(df,
        caption  = "Carbon stock (kg C/m²) by stratum and VM0033 depth interval",
        options  = list(pageLength = 20, dom = "t"),
        rownames = FALSE
      ) |>
        DT::formatRound(
          columns = intersect(
            c("mean_carbon_stock", "sd_carbon_stock", "se_carbon_stock",
              "carbon_stock_kg_m2", "mean", "sd", "se"),
            names(df)
          ),
          digits = 3
        )
    })

    output$depth_plot <- renderPlot({
      plots <- read_target("eda_plots")
      req(!is.null(plots))

      # Try common list structures returned by run_eda()
      p <- NULL
      if (is.list(plots)) {
        p <- plots$depth_profiles %||% plots[[1]]
      } else {
        p <- plots
      }
      req(!is.null(p))
      print(p)
    }, res = 96)

    output$report_links <- renderUI({
      step1_path <- file.path(project_root, "reports", "step1_nonspatial.html")
      step3_path <- file.path(project_root, "reports", "step3_random_forest.html")

      items <- list()

      if (file.exists(step1_path)) {
        items[[length(items) + 1]] <- div(class = "mb-3",
          tags$strong("Step 1 — Data quality & depth harmonization"),
          tags$br(),
          downloadButton(session$ns("dl_step1"), "Download step1_nonspatial.html",
            class = "btn btn-outline-primary btn-sm mt-1")
        )
      }

      if (file.exists(step3_path)) {
        items[[length(items) + 1]] <- div(class = "mb-3",
          tags$strong("Step 3 — Random forest prediction map"),
          tags$br(),
          downloadButton(session$ns("dl_step3"), "Download step3_random_forest.html",
            class = "btn btn-outline-primary btn-sm mt-1")
        )
      }

      if (length(items) == 0) {
        div(class = "alert alert-warning",
          "Reports are not yet available. Run Pipeline 1 first.")
      } else {
        tagList(items)
      }
    })

    output$dl_step1 <- downloadHandler(
      filename = "step1_nonspatial.html",
      content  = function(file) {
        file.copy(file.path(project_root, "reports", "step1_nonspatial.html"), file)
      }
    )

    output$dl_step3 <- downloadHandler(
      filename = "step3_random_forest.html",
      content  = function(file) {
        file.copy(file.path(project_root, "reports", "step3_random_forest.html"), file)
      }
    )
  })
}
