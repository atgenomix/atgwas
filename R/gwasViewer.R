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
  library(scales)

  ui <- fluidPage(
    titlePanel("GWAS Viewer"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Choose GWAS .assoc file",
                  accept = c(".txt", ".assoc", ".assoc.logistic", ".logistic", ".csv")),
        sliderInput("genomewideline", "Genome-wide threshold (-log10):",
                    min = 0, max = 10, value = -log10(5e-8), step = 0.1),
        sliderInput("suggestiveline", "Suggestive threshold (-log10):",
                    min = 0, max = 10, value = -log10(1e-5), step = 0.1)
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
    gwas_data <- reactive({
      req(input$file)
      df <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)
      if ("TEST" %in% colnames(df)) {
        df <- subset(df, TEST == "ADD")
      }
      if ("P" %in% colnames(df)) {
        df <- df[is.finite(df$P) & df$P > 0 & df$P <= 1, ]
      }
      df
    })

    output$table <- renderDT({
      req(gwas_data())
      gwas_data()
    }, options = list(pageLength = 20))

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

      sig   <- dat[dat$logP >= input$genomewideline, ]
      nonsig <- dat[dat$logP < input$genomewideline, ]
      chr_n <- length(unique(dat$CHR))
      colors <- hue_pal()(chr_n)

      p <- ggplot() +
        geom_point(data = nonsig,
                   aes(x = pos, y = logP, color = as.factor(CHR)),
                   size = 0.5, alpha = 0.4) +
        geom_point(data = sig,
                   aes(x = pos, y = logP, color = as.factor(CHR), text = tooltip),
                   size = 0.7) +
        scale_color_manual(values = colors) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(title = "GWAS Manhattan Plot",
             x = "Genomic Position",
             y = "-log10(P)")

      xrange <- range(dat$pos, na.rm = TRUE)
      ggplotly(p, tooltip = "text") %>%
        layout(
          hovermode = "closest",
          yaxis = list(range = c(0, max(dat$logP, na.rm = TRUE) * 1.05 + 2)),
          shapes = list(
            list(
              type  = "line",
              xref  = "x", x0 = xrange[1], x1 = xrange[2],
              yref  = "y", y0 = input$genomewideline, y1 = input$genomewideline,
              line  = list(color = "red", dash = "dash")
            ),
            list(
              type  = "line",
              xref  = "x", x0 = xrange[1], x1 = xrange[2],
              yref  = "y", y0 = input$suggestiveline, y1 = input$suggestiveline,
              line  = list(color = "blue", dash = "dot")
            )
          )
        )
    })

    observeEvent(c(input$genomewideline, input$suggestiveline), {
      dat <- base_data()[1:50000, ]
      req(dat)
      xrange <- range(dat$pos, na.rm = TRUE)
      plotlyProxy("manhattan", session) %>%
        plotlyProxyInvoke("relayout", list(
          shapes = list(
            list(
              type  = "line",
              xref  = "x", x0 = xrange[1], x1 = xrange[2],
              yref  = "y", y0 = input$genomewideline, y1 = input$genomewideline,
              line  = list(color = "red", dash = "dash")
            ),
            list(
              type  = "line",
              xref  = "x", x0 = xrange[1], x1 = xrange[2],
              yref  = "y", y0 = input$suggestiveline, y1 = input$suggestiveline,
              line  = list(color = "blue", dash = "dot")
            )
          )
        ))
    })

    output$qqplot <- renderPlot({
      dat <- base_data()[1:50000, ]
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
        pch = 19, cex = 0.5,
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
