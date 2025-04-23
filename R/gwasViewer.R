
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

    ui <- fluidPage(
      titlePanel("GWAS Viewer"),
      sidebarLayout(
        sidebarPanel(
          dbBrowserUI("dbBrowser1")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Table", shinycssloaders::withSpinner(DTOutput("table"))),
            tabPanel("QQ Plot", shinycssloaders::withSpinner(plotOutput("qqplot"))),
            tabPanel("Manhattan", shinycssloaders::withSpinner(plotOutput("manhattan")))
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

      # Reactive expression to read uploaded data
      gwas_data <- reactive({
        req(db_info$selected_db())
        print(1)
        sel_db <- db_info$selected_db()
        DBI::dbExecute(sc, paste0("USE ", sel_db))
        selected_db_name <- db_info$selected_db()
        tables <- DBI::dbListTables(sc)
        print(tables)
        req(length(tables) > 0, "No tables found in ", sel_db)
        tbl_name <- tables[1]
        tbl_name <- "output_gwas_results_delta_740715742"
        df <- DBI::dbGetQuery(sc, paste0("SELECT * FROM ", tbl_name))
        df$"P" <- as.numeric(df$"P")
        df
      })

      # Render data table
      output$table <- renderDT({
        req(gwas_data())
        gwas_data()
      }, options = list(pageLength = 20))

      # Render Manhattan plot
      output$manhattan <- renderPlot({
        dat <- gwas_data()
        req(dat$CHR, dat$BP, dat$P)
        # Convert CHR to numeric, filter invalid
        dat$CHR <- as.numeric(as.character(dat$CHR))
        dat <- dat[!is.na(dat$CHR) & is.finite(dat$BP), ]
        qqman::manhattan(
          dat,
          chr = "CHR",
          bp  = "BP",
          p   = "P",
          snp = if ("SNP" %in% names(dat)) "SNP" else NULL,
          col = c("gray30", "gray60"),
          genomewideline = -log10(5e-8),
          suggestiveline  = -log10(1e-5),
          main = "GWAS Manhattan Plot",
          ylim = c(0, max(-log10(dat$P), na.rm = TRUE) * 1.05)
        )
      })

      # Render QQ plot
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

    # Launch Shiny app
    for_run <- shinyApp(ui = ui, server = server)
    runapp <- runApp(for_run)


}

