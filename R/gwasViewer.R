
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
        sliderInput("genomewideline", "genomewide threshold:", min = 0, max = 10, value = -log10(5e-8), step = 0.1),
        sliderInput("suggestiveline", "suggestive threshold:", min = 0, max = 10, value = -log10(1e-5), step = 0.1),
        dbBrowserUI("dbBrowser1")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Table", shinycssloaders::withSpinner(DTOutput("table"))),
          tabPanel("QQ Plot", shinycssloaders::withSpinner(plotOutput("qqplot"))),
          tabPanel("Manhattan", shinycssloaders::withSpinner(plotlyOutput("manhattan")))
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

    # 快取 base plot
    manhattan_plot_base <- reactiveVal(NULL)

    observeEvent(gwas_data(), {
      dat <- gwas_data()[1:10000,]
      req(dat$CHR, dat$BP, dat$P)
      dat$CHR <- as.numeric(as.character(dat$CHR))
      dat <- dat[!is.na(dat$CHR) & is.finite(dat$BP), ]
      dat$logP <- -log10(dat$P)
      dat <- dat[order(dat$CHR, dat$BP), ]
      dat$pos <- ave(dat$BP, dat$CHR, FUN = function(x) seq_along(x))
      dat$tooltip <- paste0(
        "SNP: ", if ("SNP" %in% colnames(dat)) dat$SNP else NA,
        "<br>CHR: ", dat$CHR,
        "<br>BP: ", dat$BP,
        "<br>P: ", signif(dat$P, 3)
      )
      chr_n <- length(unique(dat$CHR))
      colors <- hue_pal()(chr_n)

      p <- ggplot(dat, aes(x = pos, y = logP, color = as.factor(CHR), text = tooltip)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = colors) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(title = "GWAS Manhattan Plot", x = "Genomic Position", y = "-log10(P)")

      manhattan_plot_base(p)
    })

    output$manhattan <- renderPlotly({
      dat <- gwas_data()[1:10000,]
      req(dat$CHR, dat$BP, dat$P)
      dat$CHR <- as.numeric(as.character(dat$CHR))
      dat <- dat[!is.na(dat$CHR) & is.finite(dat$BP), ]
      dat$logP <- -log10(dat$P)
      dat <- dat[order(dat$CHR, dat$BP), ]
      dat$pos <- ave(dat$BP, dat$CHR, FUN = function(x) seq_along(x))

      # Separate significant vs. non-significant SNPs based on genome-wide threshold
      threshold <- input$genomewideline
      sig_data    <- dat[dat$logP >= threshold, ]
      nonsig_data <- dat[dat$logP <  threshold, ]

      chr_n <- length(unique(dat$CHR))
      color_palette <- scales::hue_pal()(chr_n)

      # Build base ggplot: use non-significant points as static background
      p <- ggplot() +
        geom_point(data = nonsig_data, aes(x = pos, y = logP, color = as.factor(CHR)), size = 0.5, alpha = 0.5) +
        geom_point(data = sig_data, aes(x = pos, y = logP, color = as.factor(CHR), text = paste0(
          "CHR: ", CHR,
          "<br>BP: ", BP,
          "<br>P: ", signif(P, 3),
          if ("SNP" %in% names(dat)) paste0("<br>SNP: ", SNP) else ""
        )), size = 0.8) + # 這層給 plotly 做互動
        geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
        geom_hline(yintercept = input$suggestiveline, linetype = "dotted", color = "blue") +
        scale_color_manual(values = color_palette) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(title = "GWAS Manhattan Plot", x = "Genomic Position", y = "-log10(P)")

      # Convert only the significant points into an interactive plot
      ggplotly(p, tooltip = "text", originalData = FALSE) %>%
        layout(hovermode = "closest")
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

  for_run <- shinyApp(ui = ui, server = server)
  runapp <- runApp(for_run)
}

