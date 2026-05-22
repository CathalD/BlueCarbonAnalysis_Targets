mod_run_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "step-intro",
      h4("Run Pipeline 1 — Core Data Processing."),
      p("This step harmonizes your field cores to VM0033 standard depths, ",
        "generates exploratory plots, and (if a raster is available) fits a ",
        "random forest spatial model. Typical runtime: 5–15 minutes.")
    ),

    bslib::card(
      bslib::card_body(
        div(class = "d-flex align-items-center gap-3 mb-3",
          actionButton(ns("run_btn"), "▶ Run Pipeline 1",
            class = "btn btn-success btn-lg",
            icon  = icon("play")),
          uiOutput(ns("status_badge"))
        ),

        uiOutput(ns("progress_bar_ui")),

        div(id = ns("progress_area"),
          DT::dataTableOutput(ns("progress_table"))
        ),

        uiOutput(ns("error_detail"))
      )
    )
  )
}

mod_run_server <- function(id, project_root) {
  moduleServer(id, function(input, output, session) {

    bg_proc        <- reactiveVal(NULL)
    pipeline_status <- reactiveVal("idle")  # idle | running | done | error

    observeEvent(input$run_btn, {
      req(pipeline_status() %in% c("idle", "error"))

      pipeline_status("running")

      proc <- callr::r_bg(
        func = function(root) {
          setwd(root)
          targets::tar_make()
        },
        args       = list(root = project_root),
        supervise  = TRUE,
        stdout     = NULL,
        stderr     = NULL
      )
      bg_proc(proc)
    })

    # Poll process completion
    observe({
      req(pipeline_status() == "running")
      invalidateLater(2000, session)

      proc <- bg_proc()
      if (!is.null(proc) && !proc$is_alive()) {
        exit_code <- tryCatch(proc$get_exit_status(), error = function(e) -1L)
        pipeline_status(if (isTRUE(exit_code == 0L)) "done" else "error")
      }
    })

    # Read progress from targets store
    progress_df <- reactive({
      req(pipeline_status() != "idle")
      if (pipeline_status() == "running") invalidateLater(2000, session)

      tryCatch(
        targets::tar_progress(store = file.path(project_root, "_targets")),
        error = function(e) NULL
      )
    })

    output$status_badge <- renderUI({
      switch(pipeline_status(),
        idle    = NULL,
        running = span(class = "badge bg-primary fs-6", "⏳ Running…"),
        done    = span(class = "badge bg-success fs-6", "✓ Complete"),
        error   = span(class = "badge bg-danger  fs-6", "✗ Error — see details below")
      )
    })

    output$progress_bar_ui <- renderUI({
      df <- progress_df()
      if (is.null(df) || pipeline_status() == "idle") return(NULL)

      n_done  <- sum(df$status == "completed")
      n_total <- nrow(df)
      pct     <- if (n_total > 0) round(100 * n_done / n_total) else 0

      if (pipeline_status() == "done") pct <- 100

      bar_class <- switch(pipeline_status(),
        running = "progress-bar progress-bar-striped progress-bar-animated bg-primary",
        done    = "progress-bar bg-success",
        error   = "progress-bar bg-danger",
        "progress-bar"
      )

      div(class = "mb-3",
        div(class = "progress", style = "height: 24px;",
          div(class = bar_class,
            role  = "progressbar",
            style = paste0("width: ", pct, "%;"),
            `aria-valuenow` = pct,
            `aria-valuemin` = 0,
            `aria-valuemax` = 100,
            paste0(n_done, " / ", n_total, " targets complete")
          )
        )
      )
    })

    output$progress_table <- DT::renderDataTable({
      df <- progress_df()
      req(!is.null(df), nrow(df) > 0)

      df <- df[, intersect(c("name", "type", "status"), names(df)), drop = FALSE]

      status_icon <- function(s) {
        switch(s,
          completed  = "✓",
          started    = "⏳",
          dispatched = "⏳",
          canceled   = "—",
          errored    = "✗",
          s
        )
      }
      if ("status" %in% names(df)) {
        df$status <- paste(sapply(df$status, status_icon), df$status)
      }

      DT::datatable(df,
        options = list(
          pageLength = 20, dom = "tp",
          order = list(list(2, "desc"))
        ),
        rownames = FALSE,
        class    = "table-sm table-hover"
      )
    })

    output$error_detail <- renderUI({
      req(pipeline_status() == "error")

      err_df <- tryCatch({
        m <- targets::tar_meta(store = file.path(project_root, "_targets"))
        m[!is.na(m$error), c("name", "error"), drop = FALSE]
      }, error = function(e) NULL)

      if (is.null(err_df) || nrow(err_df) == 0) {
        return(div(class = "alert alert-danger mt-3",
          "Pipeline stopped with an error. Check the R console for details."
        ))
      }

      div(class = "alert alert-danger mt-3",
        tags$strong("Error details:"),
        tags$ul(
          lapply(seq_len(nrow(err_df)), function(i) {
            tags$li(tags$code(err_df$name[i]), ": ", err_df$error[i])
          })
        ),
        helpText("Re-running Pipeline 1 is safe — completed targets will be skipped.")
      )
    })

    reactive({
      list(
        status = pipeline_status(),
        done   = pipeline_status() == "done"
      )
    })
  })
}
