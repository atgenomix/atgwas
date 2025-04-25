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
        # fileInput("file", "Choose GWAS .assoc file",
        #   accept = c(".txt", ".assoc", ".assoc.logistic", ".logistic", ".csv")),
        sliderInput("genomewideline", "Genome-wide threshold (-log10):", min = 0, max = 10, value = -log10(5e-8), step = 0.1),
        sliderInput("suggestiveline", "Suggestive threshold (-log10):", min = 0, max = 10, value = -log10(1e-5), step = 0.1),
        dbBrowserUI("dbBrowser1")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Table", shinycssloaders::withSpinner(DTOutput("table"))),
          tabPanel("QQ Plot", shinycssloaders::withSpinner(plotOutput("qqplot"))),
          tabPanel("Manhattan", 
                   plotlyOutput("manhattan", height = "600px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    sc <- sparklyr::spark_connect(master = master, method = method, version = version)
    db_info <- dbBrowserServer("dbBrowser1", sc)

    session$onSessionEnded(function() {
      if (!is.null(sc)) {
        sparklyr::spark_disconnect(sc)
        message("Spark connection disconnected.")
      }
    })

    # gwas_data <- reactive({
    #   req(input$file)
    #   df <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)

    #   if ("TEST" %in% colnames(df)) {
    #     df <- subset(df, TEST == "ADD")
    #   }
    #   if ("P" %in% colnames(df)) {
    #     df <- df[is.finite(df$P) & df$P > 0 & df$P <= 1, ]
    #   }
    #   df[sample(nrow(df), 100000), ]  
    # })

    gwas_data <- reactive({
      req(db_info$selected_db())
      sel_db <- db_info$selected_db()
      DBI::dbExecute(sc, paste0("USE ", sel_db))
      tables <- DBI::dbListTables(sc)
      req(length(tables) > 0, "No tables found in ", sel_db)
      tbl_name <- tables[1]
      df <- DBI::dbGetQuery(sc, paste0("SELECT * FROM ", tbl_name))
      df$P <- as.numeric(df$P)
      df
    })

    output$table <- renderDT({
      req(gwas_data())
      gwas_data()
    }, options = list(pageLength = 20))


    manh_data <- reactive({
      req(gwas_data())
      prep_manhattan(gwas_data())
    })
    output$manhattan <- renderPlotly({
        info    <- manh_data()
        df2     <- info$df
        axis_df <- info$axis_df
        y_max_auto <- max(
          info$yrange[2], 
          input$genomewideline, 
          input$suggestiveline
        ) + 2

 
        chr_list <- axis_df$chr
        palette_24 <- rainbow(24)
        colors     <- palette_24[1:length(chr_list)]

        df2$color <- colors[ match(df2[[ "CHR" ]], chr_list) ]

  
        thr    <- input$genomewideline
        sig    <- df2[df2$logP >= thr, ]
        nonsig <- df2[df2$logP <  thr, ]


        fig <- plot_ly(type = "scatter", mode = "markers") %>%
          # non-significant
          add_trace(
            x         = nonsig$pos_cum,
            y         = nonsig$logP,
            marker    = list(color = nonsig$color, size = 4, opacity = 0.4),
            hoverinfo = "none",
            showlegend= FALSE
          ) %>%
          # significant
          add_trace(
            x         = sig$pos_cum,
            y         = sig$logP,
            marker    = list(color = sig$color, size = 6),
            text      = sig$tooltip,
            hoverinfo = "text",
            showlegend= FALSE
          ) %>%
          layout(
            title = "GWAS Manhattan Plot",
            xaxis = list(
              title    = "Chromosome",
              tickmode = "array",
              tickvals = axis_df$center,
              ticktext = axis_df$chr
            ),
            yaxis = list(
              title = "-log10(P)",
              range = c(0, y_max_auto)
            ),
            shapes = list(
              # genome-wide line
              list(
                type = "line", xref = "x",
                x0   = info$xrange[1], x1 = info$xrange[2],
                yref = "y",
                y0   = input$genomewideline,
                y1   = input$genomewideline,
                line = list(color = "red", dash = "dash")
              ),
              # suggestive line
              list(
                type = "line", xref = "x",
                x0   = info$xrange[1], x1 = info$xrange[2],
                yref = "y",
                y0   = input$suggestiveline,
                y1   = input$suggestiveline,
                line = list(color = "blue", dash = "dot")
              )
            ),
            hovermode = "closest"
          )

        fig
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
