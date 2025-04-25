#' @import shiny
#' @import qqman
#' @import DT
#' @import sparklyr
#' @import DBI
#' @import shinycssloaders
#'
#' @title ATGWAS
#' @description
#' Launches a Shiny application for interactive exploration of GWAS summary results
#' stored in a Spark catalog. Provides an interactive table, Manhattan plot, and QQ plot.
#' @param master String specifying the Spark master endpoint (e.g., \code{"sc://host:port"} or \code{"local"}).
#' @param method String indicating the connection method passed to \code{sparklyr::spark_connect()} (e.g., \code{"shell"} or \code{"spark_connect"}).
#' @param version String specifying the Spark version to use (default \code{"3.5"}).
#' @return Invisibly returns the Shiny application object after launching it with \code{runApp()}.
#' @export

gwasViewer <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
  library(shiny)
  library(DT)
  library(shinycssloaders)
  library(plotly)
  library(ggplot2)
  library(sparklyr)
  library(DBI)
  library(scales)

  ui <- fluidPage(
    titlePanel("GWAS Viewer"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Choose GWAS .assoc file",
          accept = c(".txt", ".assoc", ".assoc.logistic", ".logistic", ".csv")),
        sliderInput("genomewideline", "Genome-wide threshold (-log10):", min = 0, max = 10, value = -log10(5e-8), step = 0.1),
        sliderInput("suggestiveline", "Suggestive threshold (-log10):", min = 0, max = 10, value = -log10(1e-5), step = 0.1),
        # dbBrowserUI("dbBrowser1")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Table", shinycssloaders::withSpinner(DTOutput("table"))),
          tabPanel("QQ Plot", shinycssloaders::withSpinner(plotOutput("qqplot"))),
          tabPanel("Manhattan", 
                   shinycssloaders::withSpinner(plotlyOutput("manhattan", height = "600px")))
        )
      )
    )
  )

  server <- function(input, output, session) {
    # sc <- sparklyr::spark_connect(master = master, method = method, version = version)
    # db_info <- dbBrowserServer("dbBrowser1", sc)

    # session$onSessionEnded(function() {
    #   if (!is.null(sc)) {
    #     sparklyr::spark_disconnect(sc)
    #     message("Spark connection disconnected.")
    #   }
    # })

    gwas_data <- reactive({
      req(input$file)
      df <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)

      if ("TEST" %in% colnames(df)) {
        df <- subset(df, TEST == "ADD")
      }
      if ("P" %in% colnames(df)) {
        df <- df[is.finite(df$P) & df$P > 0 & df$P <= 1, ]
      }
      df[1:10000,]
    })
    # gwas_data <- reactive({
    #   req(db_info$selected_db())
    #   sel_db <- db_info$selected_db()
    #   DBI::dbExecute(sc, paste0("USE ", sel_db))
    #   tables <- DBI::dbListTables(sc)
    #   req(length(tables) > 0, "No tables found in ", sel_db)
    #   tbl_name <- tables[1]
    #   df <- DBI::dbGetQuery(sc, paste0("SELECT * FROM ", tbl_name))
    #   df$P <- as.numeric(df$P)
    #   df
    # })

    output$table <- renderDT({
      req(gwas_data())
      gwas_data()
    }, options = list(pageLength = 20))

    # Preprocess and cache interactive plot data
    base_data <- reactiveVal(NULL)
    observeEvent(gwas_data(), {
      dat <- gwas_data()
      req(dat$CHR, dat$BP, dat$P)
      dat$CHR <- as.numeric(as.character(dat$CHR))
      dat <- dat[!is.na(dat$CHR) & is.finite(dat$BP), ]
      dat$logP <- -log10(dat$P)
      dat <- dat[order(dat$CHR, dat$BP), ]
      dat$pos <- ave(dat$BP, dat$CHR, FUN = function(x) seq_along(x))
      dat$tooltip <- paste0(
        "CHR: ", dat$CHR,
        "<br>BP: ", dat$BP,
        "<br>P: ", signif(dat$P, 3),
        if ("SNP" %in% names(dat)) paste0("<br>SNP: ", dat$SNP) else ""
      )
      base_data(dat)
    })

    output$manhattan <- renderPlotly({
      dat <- base_data()
      req(dat)
      dat <- dat[1:50000, ]

      # separate sig vs non-sig
      sig    <- dat[dat$logP >= input$genomewideline, ]
      nonsig <- dat[dat$logP <  input$genomewideline, ]

      xrange <- range(dat$pos, na.rm = TRUE)
      yrange <- c(0, max(dat$logP, na.rm = TRUE) * 1.05 + 2)

      # 1) build a plotly figure from scratch
      fig <- plot_ly(
        type = "scatter", mode = "markers"
      ) %>%
        # static background points
        add_trace(
          x = nonsig$pos, y = nonsig$logP,
          marker = list(size = 4, opacity = 0.4),
          showlegend = FALSE,
          hoverinfo = "none"
        ) %>%
        # interactive significant points
        add_trace(
          x = sig$pos, y = sig$logP,
          marker = list(size = 6),
          text = sig$tooltip,
          hoverinfo = "text",
          showlegend = FALSE
        ) %>%
        layout(
          title = "GWAS Manhattan Plot",
          xaxis = list(title = "Genomic Position", range = xrange),
          yaxis = list(title = "-log10(P)", range = yrange),
          hovermode = "closest",
          shapes = list(
            # genome-wide line
            list(
              type = "line",
              xref = "x", x0 = xrange[1], x1 = xrange[2],
              yref = "y", y0 = input$genomewideline, y1 = input$genomewideline,
              line = list(color = "red", dash = "dash")
            ),
            # suggestive line
            list(
              type = "line",
              xref = "x", x0 = xrange[1], x1 = xrange[2],
              yref = "y", y0 = input$suggestiveline, y1 = input$suggestiveline,
              line = list(color = "blue", dash = "dot")
            )
          )
        )

      fig
    })

    # proxy 更新門檻線
    observeEvent(c(input$genomewideline, input$suggestiveline), {
      dat <- base_data()[1:50000, ]
      req(dat)
      xrange <- range(dat$pos, na.rm = TRUE)

      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout", list(
          shapes = list(
            list(
              type = "line",
              xref = "x", x0 = xrange[1], x1 = xrange[2],
              yref = "y", y0 = input$genomewideline, y1 = input$genomewideline,
              line = list(color = "red", dash = "dash")
            ),
            list(
              type = "line",
              xref = "x", x0 = xrange[1], x1 = xrange[2],
              yref = "y", y0 = input$suggestiveline, y1 = input$suggestiveline,
              line = list(color = "blue", dash = "dot")
            )
          )
        ))
    })


    output$qqplot <- renderPlot({
      dat <- gwas_data()
      req(dat$P)
      pvals <- dat$P
      pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
      obs <- -log10(sort(pvals))
      exp <- -log10(ppoints(length(pvals)))
      plot(
        exp, obs,
        xlab = "Expected -log10(p)",
        ylab = "Observed -log10(p)",
        main = "QQ Plot of GWAS P-values",
        pch = 19,
        cex = 0.5,
        xlim = range(exp, finite = TRUE),
        ylim = range(obs, finite = TRUE)
      )
      abline(0, 1, col = "red")
    })

    outputOptions(output, "manhattan", suspendWhenHidden = FALSE)
    outputOptions(output, "qqplot", suspendWhenHidden = FALSE)
  }

  shinyApp(ui = ui, server = server)
}
