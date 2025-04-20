
#' @import shiny
#' @import qqman
#' @import DT 


gwasViewer <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {

    ui <- fluidPage(
      titlePanel("GWAS Results Viewer"),
      sidebarLayout(
        sidebarPanel(
          fileInput("file", "Choose GWAS .assoc file",
                    accept = c(".txt", ".assoc", ".assoc.logistic", ".csv"))
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Table", DTOutput("table")),
            tabPanel("Manhattan", plotOutput("manhattan")),
            tabPanel("QQ Plot", plotOutput("qqplot"))
          )
        )
      )
    )

    server <- function(input, output, session) {
      # Reactive expression to read uploaded data
      gwas_data <- reactive({
        req(input$file)
        df <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)
        # Filter for additive model if present
        if ("TEST" %in% colnames(df)) {
          df <- subset(df, TEST == "ADD")
        }
        # Remove non-finite or invalid p-values
        if ("P" %in% colnames(df)) {
          df <- df[is.finite(df$P) & df$P > 0 & df$P <= 1, ]
        }
        df
      })

      # Render data table
      output$table <- renderDT({
        gwas_data()
      }, options = list(pageLength = 10))

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
    }

    # Launch Shiny app
    for_run <- shinyApp(ui = ui, server = server)
    runapp <- runApp(for_run)


}

