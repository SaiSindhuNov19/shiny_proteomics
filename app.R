library(shiny)
library(DT)
library(MSnbase)
library(Spectra)
library(plotly)
library(mzR)
library(MsBackendMgf)
library(MSGFplus)
library(mzID)
library(MSnID)
library(pheatmap)
library(shinyWidgets)  # For better widgets
library(shinyjs)       # For dynamic UI updates
library(tidyr)         # For spread() function
library(dplyr)         # For data manipulation

options(shiny.maxRequestSize = 100 * 1024^2)  # Increase file upload limit to 100 MB

# Define UI
ui <- navbarPage(
  "Peptide Identification Tool",
  
  # Link to external CSS file
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  
  # File Upload and Summary Page
  tabPanel(
    "File Upload",
    fluidPage(
      titlePanel("Upload mzML and FASTA Files"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          fileInput(
            "mzml_files",
            "Choose mzML Files",
            accept = c(".mzML"),
            multiple = TRUE
          ),
          fileInput(
            "fasta_file",
            "Choose FASTA File",
            accept = c(".fasta", ".fa"),
            multiple = FALSE
          ),
          tags$hr(),
          helpText("Upload mzML file and a FASTA file to perform analysis.")
        ),
        mainPanel(
          class = "main-panel",
          h4("Summary of Uploaded Files"),
          DT::dataTableOutput("fileSummary")
        )
      )
    )
  ),
  
  # Results Page
  tabPanel(
    "Visualization",
    fluidPage(
      titlePanel("Visualize mzML Files"),
      uiOutput("mzMLPlots")
    )
  ),
  
  # MSGF+ Page
  tabPanel(
    "MSGF+",
    fluidPage(
      titlePanel("Run MSGF+"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          actionButton("run_msgf", "Run MSGF+", class = "btn-primary"),
          tags$hr(),
          textOutput("msgf_output"),
          width = 3  # Adjust the width of the sidebar panel
        ),
        mainPanel(
          class = "main-panel",
          fluidRow(
            column(
              width = 12,
              h4("MSGF+ Results"),
              DT::dataTableOutput("msgf_table")
            )
          ),
          fluidRow(
            column(
              width = 12,
              h4("Interactive Heatmap of MSGF+ Results"),
              plotlyOutput("msgf_heatmap")  
            )
          ),
          width = 9  
        )
      )
    )
  ),
  
  # MSnID FDR Calculation Page
  tabPanel(
    "MSnID FDR",
    fluidPage(
      titlePanel("MSnID FDR Calculation"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          fileInput(
            "mzid_file",
            "Upload mzID File",
            accept = c(".mzid")
          ),
          numericInput(
            "msmsScore_threshold",
            "MS/MS Score Threshold (-log10(SpecEValue))",
            value = 10.0,
            min = 0
          ),
          numericInput(
            "mzError_threshold",
            "m/z Error Threshold (Thomson)",
            value = 0.1,
            min = 0
          ),
          actionButton("Filter_Peptides", "Filter Peptides", class = "btn-primary"),
          tags$hr(),
          helpText("Set thresholds and click 'filter' to filter PSMs.")
        ),
        mainPanel(
          class = "main-panel",
          h4("Filtered Peptide-Spectrum Matches (PSMs)"),
          DTOutput("filtered_table")
        )
      )
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  # Reactive expression to summarize uploaded files
  uploadedFiles <- reactive({
    mzml_files <- input$mzml_files
    fasta_file <- input$fasta_file
    
    if (is.null(mzml_files) && is.null(fasta_file)) {
      return(NULL)
    }
    
    data.frame(
      FileName = c(if (!is.null(mzml_files)) mzml_files$name else NULL,
                   if (!is.null(fasta_file)) fasta_file$name else NULL),
      FileSize = c(if (!is.null(mzml_files)) mzml_files$size / 1024 else NULL,
                   if (!is.null(fasta_file)) fasta_file$size / 1024 else NULL),
      FileType = c(rep("mzML", if (!is.null(mzml_files)) length(mzml_files$name) else 0),
                   if (!is.null(fasta_file)) "FASTA" else NULL),
      stringsAsFactors = FALSE
    )
  })
  
  # Render file summary table
  output$fileSummary <- DT::renderDataTable({
    files <- uploadedFiles()
    if (is.null(files)) {
      return(NULL)
    }
    DT::datatable(files, options = list(pageLength = 5, autoWidth = TRUE))
  })
  
  # Reactive expression for processing mzML files
  processMzML <- reactive({
    mzml_files <- input$mzml_files
    if (is.null(mzml_files)) {
      return(NULL)
    }
    
    lapply(seq_along(mzml_files$datapath), function(i) {
      file <- mzml_files$datapath[i]
      tryCatch({
        spectra_data <- Spectra(file)
        filtered_spectra <- filterIntensity(spectra_data, intensity = function(x) x >= 100)
        
        list(
          plotData = data.frame(
            mz = unlist(mz(filtered_spectra)),
            intensity = unlist(intensity(filtered_spectra))
          ),
          summary = data.frame(
            NumSpectra = length(spectra_data),
            RetentionTimeMin = round(min(rtime(spectra_data)), 2),
            RetentionTimeMax = round(max(rtime(spectra_data)), 2),
            TotalIonCurrent = round(sum(unlist(intensity(filtered_spectra))), 2),
            msLevels = paste(unique(msLevel(spectra_data)), collapse = ", ")
          )
        )
      }, error = function(e) {
        list(
          plotData = NULL,
          summary = data.frame(
            NumSpectra = NA,
            RetentionTimeMin = NA,
            RetentionTimeMax = NA,
            TotalIonCurrent = NA,
            msLevels = "Error reading file"
          )
        )
      })
    })
  })
  
  # Dynamic UI for mzML plots
  output$mzMLPlots <- renderUI({
    mzml_files <- input$mzml_files
    results <- processMzML()
    if (is.null(mzml_files) || is.null(results)) {
      return(NULL)
    }
    
    tabs <- lapply(seq_along(mzml_files$datapath), function(i) {
      file_name <- mzml_files$name[i]
      
      tabPanel(
        title = file_name,
        fluidRow(
          column(
            width = 6,
            tags$div(
              class = "plotly-chart",
              plotlyOutput(paste0("plot_", i))
            )
          ),
          column(
            width = 6,
            DT::dataTableOutput(paste0("table_", i))
          )
        ),
        fluidRow(
          column(
            width = 6,
            h4("Retention Time 3D Visualization"),
            plotOutput(paste0("rt3D_plot_", i))
          )
        )
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Render plots and tables dynamically
  observe({
    results <- processMzML()
    if (is.null(results)) {
      return()
    }
    
    for (i in seq_along(results)) {
      local({
        idx <- i
        plotData <- results[[idx]]$plotData
        summaryData <- results[[idx]]$summary
        
        output[[paste0("plot_", idx)]] <- renderPlotly({
          if (is.null(plotData)) return(NULL)
          plot_ly(
            data = plotData,
            x = ~mz,
            y = ~intensity,
            type = "scatter",
            mode = "lines"
          ) %>%
            layout(
              title = paste("Intensity Distribution Across m/z Values"),
              xaxis = list(title = "m/z"),
              yaxis = list(title = "Intensity")
            )
        })
        
        output[[paste0("table_", idx)]] <- DT::renderDataTable({
          if (is.null(summaryData)) return(NULL)
          DT::datatable(summaryData, options = list(pageLength = 5, autoWidth = TRUE))
        })
        
        # Render Retention Time 3D Visualization
        output[[paste0("rt3D_plot_", idx)]] <- renderPlot({
          mzf <- input$mzml_files$datapath[idx]
          ms <- openMSfile(mzf)
          hd <- header(ms)
          rtsel <- hd$retentionTime / 60 > 30 & hd$retentionTime / 60 < 35
          selected_scans <- which(rtsel)
          i <- selected_scans[1]
          j <- selected_scans[2]
          M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
          plot3D(M2)
        })
      })
    }
  })
  
  # MSGF+ reactive expression
  results <- reactiveValues(data = NULL)
  observeEvent(input$run_msgf, {
    req(input$mzml_files, input$fasta_file)
    tryCatch({
      # Get file paths
      mzml_files <- input$mzml_files$datapath
      fasta_file <- input$fasta_file$datapath
      file_names <- tools::file_path_sans_ext(basename(input$mzml_files$name))
      
      # Define output directory
      output_dir <- tempdir()  # Use a temporary directory
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      all_results <- list()
      all_columns <- c()
      
      # Loop through each mzML file
      for (i in seq_along(mzml_files)) {
        mzml_file <- mzml_files[i]
        file_name <- file_names[i]
        output_file <- file.path(output_dir, paste0(file_name, "_results.mzid"))
        
        # Set MSGF+ parameters
        param <- msgfPar(fasta_file)
        getMSGFpar(param)
        
        # Run MSGF+ and save output to the specified file
        runMSGF(
          param,
          mzml_file,
          msgfPath = "C:/Users/sai1.sindhu/Downloads/MSGFPlus_v20240326/MSGFPlus.jar",  
          savenames = output_file
        )
        
        # Load mzIdentML results
        ident_results <- mzID(output_file)
        
        # Extract PSM data
        psm_data <- flatten(ident_results)
        psm_data$FileName <- file_name  
        
        # Update all_columns to track unique column names
        all_columns <- union(all_columns, names(psm_data))
        
        # Add to results list
        all_results[[file_name]] <- psm_data
      }
      
      # Combine all results into a single data frame, ensuring consistent column names
      combined_results <- do.call(rbind, lapply(all_results, function(df) {
        missing_cols <- setdiff(all_columns, names(df))
        for (col in missing_cols) df[[col]] <- NA  
        return(df[, all_columns])
      }))
      
      results$data <- combined_results
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred while processing the files:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Display results in a DataTable
  output$msgf_table <- DT::renderDataTable({
    req(results$data)
    DT::datatable(
  results$data,
  extensions = 'Buttons',  # Enable buttons extension
  options = list(
    pageLength = 10,
    scrollX = TRUE,
    dom = 'Bfrtip',  # Position of buttons
    buttons = list(
      list(extend = "copy", text = "📋 Copy", className = "btn btn-info"),
      list(extend = "csv", text = "📄 CSV", className = "btn btn-success"),
      list(extend = "excel", text = "📊 Excel", className = "btn btn-primary"),
      list(extend = "pdf", text = "📕 PDF", className = "btn btn-danger"),
      list(extend = "print", text = "🖨️ Print", className = "btn btn-warning")
    )
  ),
  rownames = FALSE
)
  })
  
  # Interactive Heatmap for MSGF+ Results
  output$msgf_heatmap <- renderPlotly({
    req(results$data)
    
    # Prepare data for the heatmap
    heatmap_data <- results$data %>%
      group_by(FileName, pepseq) %>%  # Group by FileName and peptide sequence
      summarise(RawScore = sum(`ms-gf:rawscore`, na.rm = TRUE)) %>%  # Use MS-GF:RawScore
      spread(key = FileName, value = RawScore, fill = 0)  # Reshape data for heatmap
    
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_matrix) <- heatmap_data$pepseq
    
    # Create interactive heatmap using plotly
    plot_ly(
      x = colnames(heatmap_matrix),  
      y = rownames(heatmap_matrix),  
      z = heatmap_matrix,            
      type = "heatmap",
      colors = colorRamp(c("white", "blue", "red"))  
    ) %>%
      layout(
        title = "Peptide RawScore Heatmap",
        xaxis = list(title = "Samples"),
        yaxis = list(title = "Peptides")
      )
  })
  
  # MSnID FDR Calculation
  msnidObj <- reactiveVal(NULL)
  
  # Read the mzID file when uploaded
  observeEvent(input$mzid_file, {
    req(input$mzid_file)
    msnid <- MSnID(".")
    msnid <- read_mzIDs(msnid, input$mzid_file$datapath)
    msnidObj(msnid)
  })
  
  # Filter PSMs based on thresholds
  observeEvent(input$Filter_Peptides, {
    req(msnidObj())
    msnid <- msnidObj()
    
    # Assess termini and filter for fully tryptic peptides
    msnid <- assess_termini(msnid, validCleavagePattern = "[KR]\\.[^P]")
    msnid <- apply_filter(msnid, "numIrregCleavages == 0")
    
    # Add columns for filtering
    msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
    msnid$mzError <- abs(msnid$experimentalMassToCharge - msnid$calculatedMassToCharge)
    
    # Set up filter object
    filtObj <- MSnIDFilter(msnid)
    filtObj$msmsScore <- list(comparison = ">", threshold = input$msmsScore_threshold)
    filtObj$mzError <- list(comparison = "<", threshold = input$mzError_threshold)
    
    # Apply filter
    msnid <- apply_filter(msnid, filtObj)
    
    # Update the reactive value with filtered results
    msnidObj(msnid)
  })
  
  # Display filtered PSMs in a DataTable
  output$filtered_table <- DT::renderDataTable({
    req(msnidObj())
    filtered_data <- psms(msnidObj())
    
    DT::datatable(
      filtered_data, 
      options = list(
        pageLength = 10, 
        scrollX = TRUE,
        dom = 'Bfrtip',  # Enables buttons
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')  
      ), 
      rownames = FALSE, 
      extensions = 'Buttons'
    )
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
