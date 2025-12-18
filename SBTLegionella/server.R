library(shiny)
library(ggplot2)

source("global.R")
scheme <- ("/scheme/st_profiles.txt")
genes_names <- c("asd", "fla", "mip", "momp", "neu", "pil", "pro")


function(input, output, session) {

  volumes <- c(Home = fs::path_home())
  shinyDirChoose(input, "dir", roots = volumes, session = session)
  dir_path <- reactiveVal("")
  observeEvent(input$dir, {
    path <- parseDirPath(volumes, input$dir)
    if (length(path) > 0) {
      dir_path(path)
    }
  })
  output$path <- renderPrint({
    req(input$dir)
    parseDirPath(roots = volumes, input$dir)
  })

  data_base <- eventReactive(input$do_table,{
    req(input$dir)
    df <- createSampleTable(dir_path=dir_path())

    updatePickerInput(session,"sample", choices = c("All", unique(df$sample)))
    updatePickerInput(session, "gene", choices = c("All", unique(df$gene)), selected = "All")

    return (df)

  })

  filtered_data <- reactive({
    data <- data_base()
    if (input$sample != "All") {
      data <- data[data$sample %in% input$sample, ]
    }
    if (input$gene != "All") {
      data <- data[data$gene %in% input$gene, ]
    }
    return(data)
  })


  output$table <- renderDT({
    data <- data_base()
    if (input$sample != "All") {
      data <- data[data$sample %in% input$sample, ]
    }
    if (input$gene != "All") {
      data <- data[data$gene %in% input$gene, ]
    }

    datatable(
      filtered_data(),
      selection = "multiple",
      options = list(
        columnDefs = list(list(visible = FALSE, targets = c(1, 2)))
      )
    )
  })


  observeEvent(input$convert_fasta, {
    ab1tofasta(dir_path = dir_path(), table_samples = data_base())
  })

  sbt_results <- reactiveVal(NULL)

  observeEvent(input$do_sbt, {
    data <- filtered_data()
    dir_path <- dir_path()
    selected_rows <- input$table_rows_selected

    if (length(selected_rows) == 0) {
      showNotification("Select at least one row.", type = "error")
      return(NULL)
    }
    selected_data <- data[selected_rows, ]

    selected_genes <- unique(selected_data$gene)
    sbt_results <- list()


    for (i in selected_genes){
      data_gene <- selected_data[selected_data$gene == i, ]
      fastas <- data_gene[,2]
      samples <- data_gene[,4]

      db <- list.files(path = paste0(".", "/references"),
                               full.names = TRUE,
                               pattern = paste0(i,"\\.fas$"))

      if (length(db) == 0) {
        showNotification(paste0("Referencia no encontrada para el gen: ", current_gene), type = "warning")
        next
      }

      scheme <- createSchemes(dir_path(), i)
      alt <- paste0(i, "_profile(\\.txt)?$")
      profile <- list.files(path = paste0(dir_path,"/schemes"),
                            full.names = TRUE,
                            pattern = alt,
                            ignore.case = TRUE)

      if (length(profile) == 0) {
        showNotification("No profile files found in schemes folder.", type = "error")
        return(NULL)
      }

      resu <- SBT(files = fastas,
                  db = db,
                  scheme = profile)
      resu$result$sample <- samples
      sbt_results[[length(sbt_results) + 1]] <- resu$result
      print(sbt_results)

    }
    if (length(sbt_results) > 0) {
      final_results_df <- dplyr::bind_rows(sbt_results)
      if("ST" %in% colnames(final_results_df)) final_results_df$ST <- NULL
      final_results_df <- final_results_df %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(across(everything(), ~ first(na.omit(.))), .groups = 'drop')
      sbt_results(final_results_df)
    } else {
      showNotification("Ningún análisis SBT pudo completarse con éxito.", type = "error")
      sbt_results(NULL)
    }

  })

  output$sbt_table <- renderDataTable({
    req(sbt_results())
    datatable(sbt_results())
  })



  chrom_data <- eventReactive(input$do_chromatograph, {
    data <- data_base()
    req(data)
    if (input$sample != "All") {
      data <- data[data$sample == input$sample, ]
    }
    if (input$gene != "All") {
      data <- data[data$gene == input$gene, ]
    }

    selected_rows <- input$table_rows_selected
    idx <- selected_rows[1]
    file_path <- data$ab1_file[idx]
    req(!is.null(file_path), nzchar(file_path))
    return(file_path)
  })

  output$plot_chrom <- renderPlot({
    file_path <- chrom_data()
    req(file_path)

    ab <- readsangerseq(file_path)
    ab_c <- makeBaseCalls(ab, ratio = 0.33)
    chromatogram(ab_c, showcalls = "both", showhets = TRUE)
  })

}

