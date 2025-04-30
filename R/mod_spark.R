#' @title dbBrowerUI
#' @description search spark database
#' @param id module id for UI and Server
#' @return database browser UI module
#' @export

dbBrowserUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(
      inputId = ns("selected_db"),
      label   = "",
      choices = c("Loading, please wait..." = "")
    ),
  )
}

#' @title dbBrowerSever
#' @description search spark database
#' @param id module id for UI and Server
#' @param sc spark connection
#' @return database browser Server module
#' @export

dbBrowserServer <- function(id, sc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      org <- tolower(Sys.getenv("SPARK_USER"))
      print("=====org=====")
      print(org)
      print("=====org=====")
      c <- ifelse(stringr::str_equal(org, ""), "", sprintf("LIKE '*_%s'", org))
      print(c)
      db_list <- dbGetQuery(sc, sprintf("SHOW DATABASES %s", c))
      #choices <- db_list[["namespace"]]
      full_names <- db_list[["namespace"]]
      print("full_names")
      print(full_names)
      #org <- "cus_ejajocvzumxvupd"
      pattern <- paste0("_", org, "$")
      print("pattern")
      print(pattern)
      display_labels <- sub(pattern, "", full_names)
      print(display_labels)
      choices <- setNames(full_names, display_labels) 
      print(choices)
      # UI will show names, not values
      updateSelectInput(
        session,
        "selected_db",
        choices = choices,
        selected = ""
      )
    })
    selected_db <- reactive({ input$selected_db })
    return(list(
        selected_db = selected_db
    ))
  })
}




#' @title UI for connecting to spark
#' @description spark connection shiny UI
#' @param id id
#' @return pivot long format from se object
#' @export


sparkConnectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("conn_option"), "Connection Setting",
                 choices = c("Default" = "default", "Custom" = "custom"),
                 selected = "default", inline = TRUE),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'custom'", ns("conn_option")),
      textInput(ns("spark_master"), "Spark Master", value = "sc://localhost:15002"),
      textInput(ns("spark_method"), "Spark Method", value = "spark_connect"),
      textInput(ns("spark_version"), "Spark Version", value = "3.5")
    ),
    actionButton(ns("connect"), "Connect")
  )
}


#' @title Server for connecting to spark
#' @description spark connection shiny server
#' @param id id
#' @return pivot long format from se object
#' @export

sparkConnectionServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    sc <- reactiveVal(NULL)

    observeEvent(input$connect, {
      if (input$conn_option == "default") {
        master <- "local"
        print(input$conn_option)
        connection <<- sparklyr::spark_connect(master = master)
        sc(connection)
      } else {
        master <- input$spark_master
        method <- input$spark_method
        version <- input$spark_version
        connection <<- sparklyr::spark_connect(master = master, method = method, version = version )
        sc(connection)
      }
    })
    return(sc)
  })
}
