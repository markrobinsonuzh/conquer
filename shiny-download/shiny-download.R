library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## TODO: Add link to GEO/ArrayExpress data set

scrna_download_shiny <- function(data_directory, top_url) {
  ## ----------------------------------------------------------------------- ##
  ##                           Define UI                                     ##
  ## ----------------------------------------------------------------------- ##
  p_layout <- 
    shinydashboard::dashboardPage(
      skin = "blue",
      
      shinydashboard::dashboardHeader(
        title = "conquer",
        titleWidth = 200
      ),
      
      shinydashboard::dashboardSidebar(
        width = 0
      ),
        
      shinydashboard::dashboardBody(fluidRow(
        shinydashboard::tabBox(
          width = 12,
          
          tags$head(
            tags$style(HTML("
                            .shiny-output-error-validation {
                            font-size: 17px;
                            }
                            "))
            ),
          
          tabPanel("scRNA-seq data sets",
                   DT::dataTableOutput("dt_datasets"),
                   value = "select_dataset"),
          
          selected = "select_dataset"
            )
            ))
      )
  
  ## ----------------------------------------------------------------------- ##
  ##                           Define server                                 ##
  ## ----------------------------------------------------------------------- ##
  
  server_function <- function(input, output, session) {
    withProgress(message = "Retrieving data...", value = 0, {
      ## MultiAssayExperiment
      maestbl <- file.info(list.files(paste0(data_directory, "/data-mae"), pattern = "\\.rds$",
                                      full.names = TRUE))
      rownames(maestbl) <- gsub(data_directory, top_url, rownames(maestbl))
      colnames(maestbl) <- paste0("mae_", colnames(maestbl))
      maestbl$mae_link <- rownames(maestbl)
      rownames(maestbl) <- gsub("\\.rds$", "", basename(rownames(maestbl)))

      ## MultiQC
      # multiqctbl <- file.info(system(paste0("ls ", data_directory, "/data-processed/*/multiqc/*_multiqc_report.html"), 
      #                                intern = TRUE))
      multiqctbl <- file.info(list.files(paste0(data_directory, "/report-multiqc"),
                                         pattern = "_multiqc_report.html$",
                                         full.names = TRUE))
      rownames(multiqctbl) <- gsub(data_directory, top_url, rownames(multiqctbl))
      colnames(multiqctbl) <- paste0("multiqc_", colnames(multiqctbl))
      multiqctbl$multiqc_link <- rownames(multiqctbl)
      rownames(multiqctbl) <- gsub("_multiqc_report.html$", "", basename(rownames(multiqctbl)))
      
      ## scater
      scatertbl <- file.info(list.files(paste0(data_directory, "/report-scater"), 
                                        pattern = "_scater.html$",
                                        full.names = TRUE))
      rownames(scatertbl) <- gsub(data_directory, top_url, rownames(scatertbl))
      colnames(scatertbl) <- paste0("scater_", colnames(scatertbl))
      scatertbl$scater_link <- rownames(scatertbl)
      rownames(scatertbl) <- gsub("_scater.html$", "", basename(rownames(scatertbl)))

      ## Salmon
      salmontbl <- file.info(system(paste0("ls ", data_directory, "/data-processed/*/*_salmon.tar.gz"), 
                                     intern = TRUE))
      # salmontbl <- file.info(list.files(paste0(data_directory, "/data-processed"), 
      #                                   pattern = "_salmon.tar.gz$",
      #                                   full.names = TRUE, recursive = TRUE))
      rownames(salmontbl) <- gsub(data_directory, top_url, rownames(salmontbl))
      colnames(salmontbl) <- paste0("salmon_", colnames(salmontbl))
      salmontbl$salmon_link <- rownames(salmontbl)
      rownames(salmontbl) <- gsub("_salmon.tar.gz$", "", basename(rownames(salmontbl)))
      
      dtbl <- merge(merge(maestbl, multiqctbl, by = 0, all = TRUE),
                    merge(scatertbl, salmontbl, by = 0, all = TRUE),
                    by = "Row.names", all = TRUE)
      
      ## Add some basic information about the data set
      dss <- list.files(paste0(data_directory, "/data-processed"), full.names = TRUE, recursive = FALSE)
      names(dss) <- basename(dss)
      basic_info <- data.frame(t(sapply(dss, function(ds) {
        fn <- paste0(ds, "/dataset_info.txt")
        if (file.exists(fn)) {
          fn <- read.delim(fn, header = FALSE, row.names = 1, as.is = TRUE)
          c(organism = fn["organism", 1], nsamples = fn["nsamples", 1])
        } else {
          c(organism = NA, nsamples = NA)
        }
      })), stringsAsFactors = FALSE)

      dtbl <- merge(dtbl, basic_info,
                    by.x = "Row.names", by.y = 0, all = TRUE)
    })

    create_link <- function(val, title, chdate) {
      s <- sprintf('<a href="%s" target="_blank" class="btn btn-primary">Download %s</a> (%s)', 
                   val, title, chdate)
      s[is.na(val)] <- "Not available"
      s
    }
    
    output$dt_datasets <- DT::renderDataTable({
      dtbl$data_set <- dtbl$Row.names
      dtbl$MultiAssayExperiment <- create_link(dtbl$mae_link, ".rds", as.Date(dtbl$mae_mtime))
      dtbl$MultiQC_html <- create_link(dtbl$multiqc_link, ".html", as.Date(dtbl$multiqc_mtime))
      dtbl$scater_html <- create_link(dtbl$scater_link, ".html", as.Date(dtbl$scater_mtime))
      dtbl$salmon_tar <- create_link(dtbl$salmon_link, ".tar.gz", as.Date(dtbl$salmon_mtime))
      return(dtbl[, c("data_set", "organism", "nsamples", "MultiAssayExperiment", 
                      "MultiQC_html", "scater_html", "salmon_tar")])
    }, escape = FALSE)   

  }
  
  shinyApp(ui = p_layout, server = server_function)
}