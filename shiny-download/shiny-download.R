suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DT))

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
          
          tabPanel("About",
                   includeMarkdown("/home/Shared_taupo/data/seq/conquer/database/shiny-download/about_conquer.md"),
                   value = "about"),
          
          tabPanel("scRNA-seq data sets",
                   DT::dataTableOutput("dt_datasets"),
                   value = "select_dataset"),
          
          tabPanel("Changelog",
                   includeMarkdown("/home/Shared_taupo/data/seq/conquer/database/shiny-download/changelog_conquer.md"),
                   value = "changelog"),
          
          tabPanel("Excluded samples",
                   includeMarkdown("/home/Shared_taupo/data/seq/conquer/database/shiny-download/excluded_samples.md"),
                   value = "excludedsamples"),

          tabPanel("Tutorial",
                   includeMarkdown("/home/Shared_taupo/data/seq/conquer/database/shiny-download/tutorial.md"),
                   value = "tutorial"),
          
          selected = "select_dataset"
        )
      ))
    )
  
  ## ----------------------------------------------------------------------- ##
  ##                           Define server                                 ##
  ## ----------------------------------------------------------------------- ##
  
  server_function <- function(input, output, session) {
    
    withProgress(message = "Retrieving data...", value = 0, {
      
      ## ----------------- MultiAssayExperiment table ----------------------- ##
      maestbl <- file.info(list.files(paste0(data_directory, "/data-mae"), 
                                      pattern = "\\.rds$", full.names = TRUE))
      colnames(maestbl) <- paste0("mae_", colnames(maestbl))
      maestbl$mae_link <- gsub(data_directory, top_url, rownames(maestbl))
      maestbl$dataset <- gsub("(_20[0-9]{6})*\\.rds$", "", basename(rownames(maestbl)))
      maestbl <- maestbl %>% dplyr::group_by(dataset) %>% dplyr::arrange(mae_mtime) %>%
        dplyr::summarise_all(funs(last)) %>% dplyr::ungroup() %>% 
        as.data.frame()

      ## --------------------- MultiQC table -------------------------------- ##
      multiqctbl <- file.info(list.files(paste0(data_directory, "/report-multiqc"),
                                         pattern = "_multiqc_report(_20[0-9]{6})*.html$",
                                         full.names = TRUE))
      colnames(multiqctbl) <- paste0("multiqc_", colnames(multiqctbl))
      multiqctbl$multiqc_link <- gsub(data_directory, top_url, rownames(multiqctbl))
      multiqctbl$dataset <- gsub("_multiqc_report(_20[0-9]{6})*\\.html$", "", 
                                 basename(rownames(multiqctbl)))
      multiqctbl <- multiqctbl %>% dplyr::group_by(dataset) %>% 
        dplyr::arrange(multiqc_mtime) %>% dplyr::summarise_all(funs(last)) %>% 
        dplyr::ungroup() %>% as.data.frame()

      ## --------------------- scater table --------------------------------- ##
      scatertbl <- file.info(list.files(paste0(data_directory, "/report-scater"), 
                                        pattern = "_scater(_20[0-9]{6})*.html$",
                                        full.names = TRUE))
      colnames(scatertbl) <- paste0("scater_", colnames(scatertbl))
      scatertbl$scater_link <- gsub(data_directory, top_url, rownames(scatertbl))
      scatertbl$dataset <- gsub("_scater(_20[0-9]{6})*\\.html$", "", 
                                basename(rownames(scatertbl)))
      scatertbl <- scatertbl %>% dplyr::group_by(dataset) %>% 
        dplyr::arrange(scater_mtime) %>% dplyr::summarise_all(funs(last)) %>% 
        dplyr::ungroup() %>% as.data.frame()
      
      ## --------------------- Salmon table --------------------------------- ##
      salmontbl <- file.info(system(paste0("ls ", data_directory, 
                                           "/data-processed/*/*_salmon*.tar.gz"), 
                                     intern = TRUE))
      colnames(salmontbl) <- paste0("salmon_", colnames(salmontbl))
      salmontbl$salmon_link <- gsub(data_directory, top_url, rownames(salmontbl))
      salmontbl$dataset <- gsub("_salmon(_20[0-9]{6})*\\.tar.gz$", "", 
                                basename(rownames(salmontbl)))
      salmontbl <- salmontbl %>% dplyr::group_by(dataset) %>% 
        dplyr::arrange(salmon_mtime) %>% dplyr::summarise_all(funs(last)) %>% 
        dplyr::ungroup() %>% as.data.frame()
      
      ## ---------------------- Merge tables -------------------------------- ##
      dtbl <- maestbl %>% dplyr::full_join(multiqctbl, by = "dataset") %>%
        dplyr::full_join(scatertbl, by = "dataset") %>%
        dplyr::full_join(salmontbl, by = "dataset")

      ## ----------- Add some basic information about the data sets --------- ##
      dss <- list.files(paste0(data_directory, "/data-processed"), 
                        full.names = TRUE, recursive = FALSE) ## datasets
      names(dss) <- basename(dss)
      basic_info <- data.frame(t(sapply(dss, function(ds) {
        fn <- file.info(list.files(ds, pattern = "dataset_info", 
                                   full.names = TRUE)) %>%
          tibble::rownames_to_column(var = "file") %>%
          dplyr::arrange(mtime) %>% dplyr::summarise_all(funs(last)) %>% 
          as.data.frame()
        fn <- fn$file
        if (file.exists(fn)) {
          fn <- read.delim(fn, header = FALSE, row.names = 1, as.is = TRUE)
          c(organism = fn["organism", 1], ncells = fn["nsamples", 1],
            pmid = fn["PMID", 1], datalink = fn["datalink", 1], 
            ID = fn["shortname", 1], description = fn["description", 1])
        } else {
          c(organism = NA, ncells = NA, pmid = NA, datalink = NA, ID = NA,
            description = NA)
        }
      })), stringsAsFactors = FALSE) %>% 
        tibble::rownames_to_column(var = "dataset")

      dtbl <- dtbl %>% dplyr::full_join(basic_info, by = "dataset")
    })

    ## ----------------- Help functions to create links --------------------- ##
    create_link <- function(val, title, chdate) {
      s <- sprintf('<a href="%s" target="_blank" class="btn btn-primary">Download %s</a> (%s)', 
                   val, title, chdate)
      s[is.na(val)] <- "Not available"
      s
    }
    
    create_link_2 <- function(val, title, val2) {
      s1 <- sapply(seq_along(val), function(i) {
        if (!is.na(val[i])) sprintf('<a href="%s" target="_blank"> %s</a>',
                                    val[i], title[i])
        else title[i]
      })
      s2 <- sapply(seq_along(val), function(i) {
        if (!is.na(val2[i])) sprintf('<a href="%s" target="_blank"> %s</a>',
                                     paste0("http://www.ncbi.nlm.nih.gov/pubmed/", val2[i]),
                                     paste0(" (PMID ", val2[i], ")"))
        else ""
      })
      paste0(s1, s2)
    }
    
    ## --------------------- Final data table ------------------------------- ##
    output$dt_datasets <- DT::renderDataTable({
      dtbl$`Data set` <- create_link_2(dtbl$datalink, dtbl$dataset, dtbl$pmid)
      dtbl$MultiAssayExperiment <- create_link(dtbl$mae_link, ".rds", as.Date(dtbl$mae_mtime))
      dtbl$`MultiQC report` <- create_link(dtbl$multiqc_link, ".html", as.Date(dtbl$multiqc_mtime))
      dtbl$`scater report` <- create_link(dtbl$scater_link, ".html", as.Date(dtbl$scater_mtime))
      dtbl$`salmon archive` <- create_link(dtbl$salmon_link, ".tar.gz", as.Date(dtbl$salmon_mtime))
      dtbl$`Number of cells` <- as.numeric(as.character(dtbl$ncells))
      dtbl$`Brief description` <- dtbl$description
      # return(dtbl[, c("Data set", "ID", "organism", "ncells", "MultiAssayExperiment", 
      #                 "MultiQC report", "scater report", "salmon archive")])
      return(datatable(dtbl[, c("Data set", "ID", "organism", "Brief description",
                                "Number of cells", "MultiAssayExperiment", 
                                "MultiQC report", "scater report", "salmon archive")],
                       escape = FALSE, 
                       options = list(scrollX = TRUE)))
    }, escape = FALSE)   

  }
  
  shinyApp(ui = p_layout, server = server_function)
}
