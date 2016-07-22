library(MultiAssayExperiment)
library(Biostrings)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(tximport)
library(GEOquery)
library(rmarkdown)
library(Rtsne)

generate_report <- function(id, maex, phenoid, output_format = NULL,
                            output_file = NULL, output_dir = "./", ...){
  ## This function was written by Nicholas Hamilton and obtained from 
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r
  
  ## Give the path to the template file
  theFile <- "scater_template.Rmd"
  
  ## Process the arguments
  args <- list(...)
  args$input <- theFile
  args$output_dir <- output_dir
  args$output_format <- output_format
  args$output_file <- output_file
  
  ## Render the report
  outputFileName <- do.call('render', args = args)
  invisible(outputFileName)
}

#' Function to process fastq files
#' 
#' fastq files will be streamed into FastQC
#' (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and Salmon
#' (https://combine-lab.github.io/salmon/) for analysis. Results will be
#' summarized using MultiQC (http://multiqc.info/) and scater
#' (https://www.bioconductor.org/packages/release/bioc/html/scater.html).
#' 
#' @param id A data set ID (typically a GSE ID). It is assumed that there is a 
#'   subfolder of \code{./data-raw} with this name, and that it contains a file 
#'   named \code{{ID}_SraRunInfo.csv}.
#' @param rtype Either "single" or "paired"
#' @param index The path to the Salmon index
#' @param libtype The \code{LIBTYPE} argument passed to Salmon
#' @param salmonbin The path to the Salmon binary
#' @param fastqcbin The path to the FastQC binary
#' @param sncol The column in the SraRunInfo file corresponding to sample names
#' @param txgenemap A data.frame mapping between transcripts and genes
#' @param geodata TRUE or FALSE, whether the data set is from GEO or not. If
#'   TRUE, the function assumes that there is a file
#'   \code{./data-raw/{ID}_series_matrix.txt.gz}. If FALSE, you have to provide
#'   a \code{phenofile}.
#' @param If the data set is not a GEO data set, a text file with phenotype
#'   information
#' @param gene_granges A GRanges object with gene information
#' @param tx_granges A GRanges object with transcript information
#' @param groupid The name of the annotation column (from the phenotype data
#'   file) that will be used to color and stratify samples in the scater
#'   analysis
#' @param organism The organism
#' @param genome The genome build
#' 
#' @return Does not return anything, but saves processed files and reports in
#'   the ./data-processed and ./qc-scater subfolders
#'   
#' @author Charlotte Soneson
#' 
process_data <- function(id, rtype, index, libtype, salmonbin = "salmon",
                         fastqcbin = "fastqc", sncol = "SampleName", 
                         txgenemap = NULL, geodata = TRUE,
                         phenofile = NULL, gene_granges = NULL, 
                         tx_granges = NULL, groupid = NULL,
                         organism, genome) {

  ## Read run info downloaded from SRA
  x <- read.delim(paste0("data-raw/", id, "/", id, "_SraRunInfo.csv"), 
                  header = TRUE, as.is = TRUE, sep = ",")
  samples <- unique(x[, sncol])

  ## ------------------------------------------------------------------------ ##
  ##                              FastQC                                      ##
  ## ------------------------------------------------------------------------ ##
  if (rtype == "single") {
    for (smp in samples) {
      ## Find all the SRA runs corresponding to this sample. They will be merged
      ## together when running FastQC
      runs <- x$Run[x[, sncol] == smp]
      files <- paste(paste0("<(./stream_ena ", runs, ".fastq)"), collapse = " ")
      
      if (!file.exists(paste0("data-processed/", id, "/fastqc/", smp))) {
        mkd <- sprintf("mkdir %s",
                       paste0("data-processed/", id, "/fastqc/", smp))
        system(mkd)
        
        fastqc <- sprintf("bash -c 'cat %s | %s --noextract -o %s -f fastq stdin:%s'",
                          files,
                          fastqcbin, 
                          paste0("data-processed/", id, "/fastqc/", smp),
                          smp)
        system(fastqc)
      } else {
        message("FastQC has already been run for ", smp)
      }
    }
  } else if (rtype == "paired") {
    for (smp in samples) {
      ## Find all the SRA runs corresponding to this sample. They will be merged
      ## together when running FastQC
      runs <- x$Run[x[, sncol] == smp]
      files1 <- paste(paste0("<(./stream_ena ", runs, "_1.fastq)"), collapse = " ")
      files2 <- paste(paste0("<(./stream_ena ", runs, "_2.fastq)"), collapse = " ")
      
      if (!file.exists(paste0("data-processed/", id, "/fastqc/", smp))) {
        mkd <- sprintf("mkdir -p %s",
                       paste0("data-processed/", id, "/fastqc/", smp))
        system(mkd)
        
        ## Read 1
        fastqc <- sprintf("bash -c 'cat %s | %s --noextract -o %s -f fastq stdin:%s'",
                          files1,
                          fastqcbin, 
                          paste0("data-processed/", id, "/fastqc/", smp),
                          paste0(smp, "_1"))
        system(fastqc)

        ## Read 2
        fastqc <- sprintf("bash -c 'cat %s | %s --noextract -o %s -f fastq stdin:%s'",
                          files2,
                          fastqcbin, 
                          paste0("data-processed/", id, "/fastqc/", smp),
                          paste0(smp, "_2"))
        system(fastqc)
      } else {
        message("FastQC has already been run for ", smp)
      }
    }
  }

  ## ------------------------------------------------------------------------ ##
  ##                              Salmon                                      ##
  ## ------------------------------------------------------------------------ ##
  if (rtype == "single") {
    for (smp in samples) {
      ## Find all the SRA runs corresponding to this sample. They will be merged
      ## together when running Salmon
      runs <- x$Run[x[, sncol] == smp]
      files <- paste(paste0("<(./stream_ena ", runs, ".fastq)"), collapse = " ")
      
      if (!file.exists(paste0("data-processed/", id, "/salmon/", smp))) {
        salmon <- sprintf("bash -c '%s quant -l %s -i %s -r <(cat %s) -o %s'",
                          salmonbin, 
                          libtype,
                          index,
                          files, 
                          paste0("data-processed/", id, "/salmon/", smp))
        system(salmon)
      } else {
        message("Salmon has already been run for ", smp)
      }
    }
  } else if (rtype == "paired") {
    for (smp in samples) {
      ## Find all the SRA runs corresponding to this sample. They will be merged
      ## together when running Salmon
      runs <- x$Run[x[, sncol] == smp]
      files1 <- paste(paste0("<(./stream_ena ", runs, "_1.fastq)"), collapse = " ")
      files2 <- paste(paste0("<(./stream_ena ", runs, "_2.fastq)"), collapse = " ")
      
      if (!file.exists(paste0("data-processed/", id, "/salmon/", smp))) {
        salmon <- sprintf("bash -c '%s quant -l %s -i %s -1 <(cat %s) -2 <(cat %s) -o %s'",
                          salmonbin, 
                          libtype,
                          index,
                          files1,
                          files2, 
                          paste0("data-processed/", id, "/salmon/", smp))
        system(salmon)
      } else {
        message("Salmon has already been run for ", smp)
      }
    }
  }
  ## Compress all Salmon output in a tar archive
  targz <- sprintf("bash -c 'tar -C data-processed/ -czf %s %s'",
                   paste0("data-processed/", id, "/", id, "_salmon.tar.gz"),
                   paste0(id, "/salmon"))
  system(targz)

  ## ------------------------------------------------------------------------ ##
  ##                             MultiQC                                      ##
  ## ------------------------------------------------------------------------ ##
  mqc <- sprintf("bash -c 'multiqc -o %s -n %s -f %s'",
                 paste0("report-multiqc"),
                 paste0(id, "_multiqc_report.html"),
                 paste0("data-processed/", id))
  system(mqc)
  
  ## ------------------------------------------------------------------------ ##
  ##                             tximport                                     ##
  ## ------------------------------------------------------------------------ ##
  files <- paste0(list.files(paste0("data-processed/", id, "/salmon"), 
                             full.names = TRUE), "/quant.sf")
  names(files) <- basename(gsub("/quant.sf", "", files))
  txi_tx <- tximport(files = files, type = "salmon", txIn = TRUE, txOut = TRUE)
  txi_gene <- summarizeToGene(txi = txi_tx, tx2gene = txgenemap,
                              countsFromAbundance = "no")
  txi_gene_lstpm <- summarizeToGene(txi = txi_tx, tx2gene = txgenemap, 
                                    countsFromAbundance = "lengthScaledTPM")

  if (geodata == TRUE) {
    geo <- getGEO(filename = paste0("data-raw/", id, "/", id, "_series_matrix.txt.gz"),
                  getGPL = FALSE)
    meta <- pData(geo)
  } else {
    meta <- read.delim(phenofile, header = TRUE, row.names = 1, as.is = TRUE)
  }
  stopifnot(all(colnames(txi_tx$counts) %in% rownames(meta)))
  meta <- meta[match(colnames(txi_tx$counts), rownames(meta)), ]
  
  ## ------------------------------------------------------------------------ ##
  ##                  Generate MultiAssayExperiment                           ##
  ## ------------------------------------------------------------------------ ##
  stopifnot(all(rownames(txi_gene_lstpm$counts) == rownames(txi_gene$counts)))
  stopifnot(all(colnames(txi_gene_lstpm$counts) == colnames(txi_gene$counts)))
  
  generse <- SummarizedExperiment(assays = list(TPM = txi_gene$abundance,
                                                count = txi_gene$counts,
                                                count_lstpm = txi_gene_lstpm$counts,
                                                avetxlength = txi_gene$length),
                                  rowRanges = gene_granges[rownames(txi_gene$abundance)])
  
  txrse <- SummarizedExperiment(assays = list(TPM = txi_tx$abundance,
                                              count = txi_tx$counts,
                                              efflength = txi_tx$length),
                                rowRanges = tx_granges[rownames(txi_tx$abundance)])
  
  ## Generate MultiAssayExperiment
  mae <- MultiAssayExperiment(ExperimentList = list(gene = generse,
                                                    tx = txrse),
                              pData = droplevels(meta))
  mae@metadata <- list(genome = genome, 
                       organism = organism,
                       index = basename(index))
  
  saveRDS(mae, file = paste0("data-mae/", id, ".rds"))
  
  ## ------------------------------------------------------------------------ ##
  ##                  Write basic information to file                         ##
  ## ------------------------------------------------------------------------ ##
  infodf <- data.frame(nsamples = nrow(pData(mae)),
                       organism = organism,
                       genome = genome,
                       ntranscripts = nrow(txi_tx$counts),
                       ngenes = nrow(txi_gene$counts))
  write.table(t(infodf), file = paste0("data-processed/", id, "dataset_info.txt"),
              row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  ## ------------------------------------------------------------------------ ##
  ##                    Generate scater QC report                             ##
  ## ------------------------------------------------------------------------ ##
  generate_report(id = id, maex = mae, phenoid = groupid, 
                  output_format = "html_document",
                  output_file = paste0(id, "_scater.html"),
                  output_dir = paste0("report-scater"))
}

