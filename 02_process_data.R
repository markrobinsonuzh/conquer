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
library(rjson)
source("00_help_functions.R")

generate_report <- function(id, maex, phenoid, output_format = NULL,
                            output_file = NULL, output_dir = "./", nrw, lps, ...){
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
#' @param phenofile If the data set is not a GEO data set, a text file with
#'   phenotype information
#' @param gene_granges A GRanges object with gene information
#' @param tx_granges A GRanges object with transcript information
#' @param groupid The name of the annotation column (from the phenotype data 
#'   file) that will be used to color and stratify samples in the scater 
#'   analysis
#' @param organism The organism
#' @param genome The genome build
#' @param dotrim Whether or not to also run analysis on trimmed data (using 
#'   cutadapt)
#' @param adapterseq The sequence that will be provided to cutadapt (via the 
#'   \code{-a} argument)
#' @param cutadaptbin The path to the cutadapt binary
#' @param tmp_dir The temporary folder where fastq files will be stored during
#'   processing
#' @param nrw The number of rows to split the figure legends over
#' @param lps The position of the legend
#' @param pmid A PubMed ID that can be linked to the dataset
#' @param datalink A URL where the dataset can be found
#' @param shortname An informative identifier for the dataset
#' @param multiqcbin The path to the multiqc binary
#' @param aspects Which parts of the data processing that will be run
#' @param topdir The top directory for reading and writing data. It is assumed
#'   that this directory has subdirectories data-mae, data-processed, data-raw,
#'   report-multiqc and report-scater
#'   
#' @return Does not return anything, but saves processed files and reports in 
#'   the ./data-processed, ./data-mae, ./report-multiqc and ./report-scater
#'   subfolders
#'   
#' @author Charlotte Soneson
#'   
process_data <- function(id, rtype, index, libtype, salmonbin = "salmon",
                         fastqcbin = "fastqc", sncol = "SampleName", 
                         txgenemap = NULL, geodata = TRUE,
                         phenofile = NULL, gene_granges = NULL, 
                         tx_granges = NULL, groupid = NULL,
                         organism, genome, dotrim = FALSE, adapterseq = NULL,
                         cutadaptbin = "cutadapt", tmp_dir = "tmp", nrw = NULL,
                         lps = "right", pmid = NA, datalink = NA, shortname = NA, 
                         multiqcbin = multiqcbin, 
                         aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"),
                         topdir = ".") {

  ## Read run info downloaded from SRA
  x <- read.delim(paste0(topdir, "/data-raw/", id, "/", id, "_SraRunInfo.csv"), 
                  header = TRUE, as.is = TRUE, sep = ",")
  samples <- unique(x[, sncol])

  if (any(c("fastqc", "salmon") %in% aspects)) {
    for (smp in samples) {
      ## Find all the SRA runs corresponding to this sample. They will be merged
      ## together in the analysis
      runs <- x$Run[x[, sncol] == smp]
      
      ## Put together a file list
      if (rtype == "single") {
        files <- paste(paste0("<(./stream_ena ", runs, ".fastq)"), collapse = " ")
      } else if (rtype == "paired") {
        files1 <- paste(paste0("<(./stream_ena ", runs, "_1.fastq)"), collapse = " ")
        files2 <- paste(paste0("<(./stream_ena ", runs, "_2.fastq)"), collapse = " ")
      } 
      
      if (rtype == "single") {
        if (!all(file.exists(paste0(topdir, "/data-processed/", id, "/fastqc/", smp, "/", smp, "_fastqc.html"),
                             paste0(topdir, "/data-processed/", id, "/salmon/", smp, "/quant.sf")))) {
          
          ## Download fastq file and save temporarily
          message("Downloading fastq file for ", smp)
          dwl <- sprintf("bash -c 'cat %s > %s'",
                         files,
                         paste0(tmp_dir, "/", smp, ".fastq"))
          system(dwl)
          ## Trim 
          if (dotrim) {
            message("Trimming fastq file for ", smp)
            trim_single(cutadapt_dir = paste0(topdir, "/data-processed/", id, "/cutadapt"), 
                        smp = smp, adapterseq = adapterseq, 
                        cutadaptbin = cutadaptbin, tmp_dir = tmp_dir)
            files <- paste0(tmp_dir, "/", smp, ".trim.fastq")
          } else {
            files <- paste0(tmp_dir, "/", smp, ".fastq")
          }
          if ("fastqc" %in% aspects)
            fastqc_single(fastqc_dir = paste0(topdir, "/data-processed/", id, "/fastqc"), 
                          smp = smp, files = files, fastqcbin = fastqcbin, appd = "") 
          if ("salmon" %in% aspects)
            salmon_single(salmon_dir = paste0(topdir, "/data-processed/", id, "/salmon"), 
                          smp = smp, files = files, salmonbin = salmonbin, 
                          libtype = libtype, index = index)
        } else {
          message("Output files for ", smp, " already exist.")
        }
        if (file.exists(files)) unlink(files)
      } else if (rtype == "paired") {
        if (!all(file.exists(paste0(topdir, "/data-processed/", id, "/fastqc/",
                                    smp, "/", smp, c("_1", "_2"), "_fastqc.html"),
                             paste0(topdir, "/data-processed/", id, "/salmon/", 
                                    smp, "/quant.sf")))) {
          ## Download fastq files and save temporarily
          message("Downloading fastq files for ", smp)
          dwl1 <- sprintf("bash -c 'cat %s > %s'",
                          files1,
                          paste0(tmp_dir, "/", smp, "_1.fastq"))
          system(dwl1)
          dwl2 <- sprintf("bash -c 'cat %s > %s'",
                          files2,
                          paste0(tmp_dir, "/", smp, "_2.fastq"))
          system(dwl2)
          ## Trim
          if (dotrim) {
            message("Trimming fastq files for ", smp)
            trim_paired(cutadapt_dir = paste0(topdir, "/data-processed/", id, "/cutadapt"), 
                        smp = smp, adapterseq = adapterseq, cutadaptbin = cutadaptbin,
                        tmp_dir = tmp_dir)
            files1 <- paste0(tmp_dir, "/", smp, "_1.trim.fastq")
            files2 <- paste0(tmp_dir, "/", smp, "_2.trim.fastq")
          } else {
            files1 <- paste0(tmp_dir, "/", smp, "_1.fastq")
            files2 <- paste0(tmp_dir, "/", smp, "_2.fastq")
          }
          if ("fastqc" %in% aspects)
            fastqc_paired(fastqc_dir = paste0(topdir, "/data-processed/", id, "/fastqc"), 
                          smp = smp, files1 = files1, files2 = files2, 
                          fastqcbin = fastqcbin)
          if ("salmon" %in% aspects)
            salmon_paired(salmon_dir = paste0(topdir, "/data-processed/", id, "/salmon"), 
                          smp = smp, files1 = files1, files2 = files2, 
                          salmonbin = salmonbin, libtype = libtype, index = index)
        } else {
          message("Output files for ", smp, " already exist.")
        }
        if (file.exists(files1)) unlink(files1)
        if (file.exists(files2)) unlink(files2)
      }
    }
  }
  
  if ("salmon" %in% aspects) {
    ## Compress all Salmon output in a tar archive
    message("Compressing Salmon output for ", id)
    targz <- sprintf("bash -c 'tar -C %s/data-processed/ -czf %s %s'",
                     topdir, 
                     paste0(topdir, "/data-processed/", id, "/", id, "_salmon.tar.gz"),
                     paste0(id, "/salmon"))
    system(targz)
    
    ## Create summary table from Salmon parameters and results
    message("Creating Salmon summary table for ", id)
    smps <- list.files(paste0(topdir, "/data-processed/", id, "/salmon"), full.names = TRUE)
    names(smps) <- basename(smps)
    summary_table_salmon <- as.data.frame(t(sapply(smps, function(s) {
      cmdinfo <- fromJSON(file = paste0(s, "/cmd_info.json"))
      cmdinfo <- c(salmon_version = cmdinfo[["salmon_version"]],
                   libtype = cmdinfo[["libType"]],
                   index = basename(cmdinfo[["index"]]))
      if (file.exists(paste0(s, "/aux_info"))) {
        metainfo <- fromJSON(file = paste0(s, "/aux_info/meta_info.json"))
      } else {
        metainfo <- fromJSON(file = paste0(s, "/aux/meta_info.json"))
      }
      metainfo <- c(num_processed = metainfo[["num_processed"]],
                    num_mapped = metainfo[["num_mapped"]],
                    percent_mapped = round(metainfo[["percent_mapped"]], 3))
      c(cmdinfo, metainfo)
    })))
    write.table(cbind(sample = rownames(summary_table_salmon), summary_table_salmon), 
                file = paste0(topdir, "/data-processed/", id, "/summary_table_salmon.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }

  ## ------------------------------------------------------------------------ ##
  ##                             MultiQC                                      ##
  ## ------------------------------------------------------------------------ ##
  if ("multiqc" %in% aspects) {
    message("Running MultiQC for ", id)
    mqc <- sprintf("bash -c '%s -o %s -n %s -f %s'",
                   multiqcbin, 
                   paste0(topdir, "/report-multiqc"),
                   paste0(id, "_multiqc_report.html"),
                   paste0(topdir, "/data-processed/", id))
    system(mqc)
  }
  
  if ("mae" %in% aspects) {
    ## ------------------------------------------------------------------------ ##
    ##                             tximport                                     ##
    ## ------------------------------------------------------------------------ ##
    message("Reading expression levels for ", id)
    files <- paste0(list.files(paste0(topdir, "/data-processed/", id, "/salmon"), 
                               full.names = TRUE), "/quant.sf")
    names(files) <- basename(gsub("/quant.sf", "", files))
    txi_tx <- tximport(files = files, type = "salmon", txIn = TRUE, txOut = TRUE)
    txi_gene <- summarizeToGene(txi = txi_tx, tx2gene = txgenemap,
                                countsFromAbundance = "no")
    txi_gene_lstpm <- summarizeToGene(txi = txi_tx, tx2gene = txgenemap, 
                                      countsFromAbundance = "lengthScaledTPM")
    
    if (geodata == TRUE) {
      geo <- getGEO(filename = paste0(topdir, "/data-raw/", id, "/", id, "_series_matrix.txt.gz"),
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
    mae <- MultiAssayExperiment(experiments = list(gene = generse,
                                                   tx = txrse),
                                pData = droplevels(meta))
    mae@metadata <- list(genome = genome, 
                         organism = organism,
                         index = basename(index))
    
    saveRDS(mae, file = paste0(topdir, "/data-mae/", id, ".rds"))
    
    ## ------------------------------------------------------------------------ ##
    ##                  Write basic information to file                         ##
    ## ------------------------------------------------------------------------ ##
    infodf <- data.frame(nsamples = nrow(pData(mae)),
                         organism = organism,
                         genome = genome,
                         ntranscripts = nrow(txi_tx$counts),
                         ngenes = nrow(txi_gene$counts),
                         PMID = pmid,
                         datalink = datalink,
                         shortname = shortname)
    write.table(t(infodf), file = paste0(topdir, "/data-processed/", id, "/dataset_info.txt"),
                row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
  
  ## ------------------------------------------------------------------------ ##
  ##                    Generate scater QC report                             ##
  ## ------------------------------------------------------------------------ ##
  if ("scater" %in% aspects) {
    mae <- readRDS(paste0(topdir, "/data-mae/", id, ".rds"))
    message("Generating scater report for ", id)
    generate_report(id = id, maex = mae, phenoid = groupid, 
                    output_format = "html_document",
                    output_file = paste0(id, "_scater.html"),
                    output_dir = paste0(topdir, "/report-scater"),
                    nrw = nrw, lps = lps)
  }
}

