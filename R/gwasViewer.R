#' @import shiny
#' @import qqman
#' @import DT
#' @import sparklyr
#' @import DBI
#' @import shinycssloaders
#' @import manhattanly
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
  # 如果沒裝過，先：
  # if (!requireNamespace("BiocManager", quietly=TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install("manhattanly")

  library(shiny)
  library(DT)
  library(shinycssloaders)
  library(plotly)
  library(manhattanly)

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
          tabPanel("Manhattan", plotlyOutput("manhattan", height = "600px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    # 1. 讀資料並簡單過濾
    gwas_data <- reactive({
      req(input$file)
      df <- read.table(input$file$datapath,
                      header = TRUE, stringsAsFactors = FALSE)
      if ("TEST" %in% colnames(df)) df <- subset(df, TEST == "ADD")
      if ("P"    %in% colnames(df)) df <- df[is.finite(df$P) & df$P > 0 & df$P <= 1, ]
      df
    })

    # 2. 資料表
    output$table <- renderDT({
      req(gwas_data())
      gwas_data()
    }, options = list(pageLength = 20))

    # 3. QQ Plot（保持不變）
    output$qqplot <- renderPlot({
      dat <- gwas_data()
      req(dat$P)
      pvals <- dat$P
      pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
      obs <- -log10(sort(pvals))
      exp <- -log10(ppoints(length(pvals)))
      plot(exp, obs,
          xlab = "Expected -log10(P)",
          ylab = "Observed -log10(P)",
          main = "QQ Plot of GWAS P-values",
          pch = 19, cex = 0.5,
          xlim = range(exp, finite = TRUE),
          ylim = range(obs, finite = TRUE))
      abline(0, 1, col = "red")
    })

    # 4. 用 manhattanly 畫 Interactive Manhattan
    output$manhattan <- renderPlotly({
      dat <- gwas_data()
      req(dat)
      manhattanly::manhattanly(
        dat,
        chr            = "CHR",
        bp             = "BP",
        p              = "P",
        snp            = if ("SNP" %in% names(dat)) "SNP" else NULL,
        col            = c("gray30", "gray60"),
        genomewideline = input$genomewideline,
        suggestiveline = input$suggestiveline,
        title          = "Interactive GWAS Manhattan"
      )
    })

    # 確保即使面板隱藏也持續運作
    outputOptions(output, "manhattan", suspendWhenHidden = FALSE)
    outputOptions(output, "qqplot",    suspendWhenHidden = FALSE)
  }

  shinyApp(ui = ui, server = server)
}
