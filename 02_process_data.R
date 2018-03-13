suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(countsimQC))
source("00_help_functions.R")
source("05_umi_functions.R")

#' Function to process fastq files
#' 
#' fastq files will be streamed into FastQC 
#' (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and Salmon 
#' (https://combine-lab.github.io/salmon/) or RapMap
#' (https://github.com/COMBINE-lab/RapMap) + umis(https://github.com/vals/umis)
#' for analysis. Results will be summarized using MultiQC (http://multiqc.info/)
#' and scater 
#' (https://www.bioconductor.org/packages/release/bioc/html/scater.html).
#' 
#' @param id A data set ID (typically a GSE ID). It is assumed that there is a 
#'   subfolder of \code{./data-raw} with this name, and that it contains a file 
#'   named \code{{ID}_SraRunInfo.csv}.
#' @param dtype Either "fl" (full length) or "umi"
#' @param rtype Either "single" or "paired"
#' @param organism The organism
#' @param genome The genome build
#' @param pmid A PubMed ID that can be linked to the dataset
#' @param datalink A URL where the dataset can be found
#' @param shortname An informative identifier for the dataset
#' @param description A very brief description of the data set to go in the
#'   repository table
#' @param protocol The protocol that was used to process the cells
#' @param protocoltype Either "full-length" or "UMI", the type of protocol used
#'   to process the cells
#' @param dotrim Whether or not to run analysis on trimmed data (using 
#'   cutadapt)
#' @param cutadaptbin The path to the cutadapt binary
#' @param adapterseq The sequence that will be provided to cutadapt (via the 
#'   \code{-a} argument)
#' @param fastqcbin The path to the FastQC binary
#' @param multiqcbin The path to the multiqc binary
#' @param salmonbin The path to the Salmon binary
#' @param salmonindex The path to the Salmon index
#' @param kallistobin The path to the Salmon binary
#' @param kallistoindex The path to the Salmon index
#' @param libtype The \code{LIBTYPE} argument passed to Salmon
#' @param bias Whether to run Salmon with bias correction
#' @param rapmapbin The path to the RapMap binary
#' @param rapmapindex The path to the RapMap index
#' @param umis_transform The path to the transform file (in JSON format) for tag
#'   counting with umis
#' @param sncol The column in the SraRunInfo file corresponding to sample names
#' @param groupid The name of the annotation column (from the phenotype data 
#'   file) that will be used to color and stratify samples in the scater 
#'   analysis
#' @param geodata TRUE or FALSE, whether the data set is from GEO or not. If 
#'   TRUE, the function assumes that there is a file 
#'   \code{./data-raw/{ID}_series_matrix.txt.gz}. If FALSE, you have to provide 
#'   a \code{phenofile}.
#' @param phenofile If the data set is not a GEO data set, a text file with
#'   phenotype information
#' @param gene_granges A GRanges object with gene information
#' @param tx_granges A GRanges object with transcript information
#' @param txgenemap A data.frame mapping between transcripts and genes
#' @param tmpdir The temporary folder where fastq files will be stored during
#'   processing
#' @param topdir The top directory for reading and writing data. It is assumed
#'   that this directory has subdirectories data-mae, data-processed, data-raw,
#'   report-multiqc and report-scater
#' @param nrw The number of rows to split the figure legends over
#' @param lps The position of the legend
#' @param aspects Which parts of the data processing that will be run
#' @param force Whether to force recalculation of FastQC/Salmon results even in
#'   cases where the files already exist.
#'   
#' @return Does not return anything, but saves processed files and reports in 
#'   the data-processed, data-mae, report-multiqc and report-scater subfolders
#'   of topdir.
#'   
#' @author Charlotte Soneson
#'   
process_data <- function(id, dtype, rtype, organism, genome, 
                         pmid = NA, datalink = NA, shortname = NA, 
                         description = "", protocol = "", protocoltype = "",
                         dotrim = FALSE, cutadaptbin, adapterseq = NULL,
                         fastqcbin, multiqcbin, 
                         salmonbin, salmonindex, 
                         kallistobin, kallistoindex,
                         libtype, bias = FALSE, 
                         rapmapbin, rapmapindex, umis_transform, cell_barcodes, 
                         sncol = "SampleName", groupid,
                         geodata = TRUE, phenofile = NULL, 
                         gene_granges = NULL, tx_granges = NULL, txgenemap = NULL, 
                         tmpdir = "tmp", topdir = ".", 
                         nrw = NULL, lps = "right", 
                         aspects = c("fastqc", "salmon", "multiqc", "mae", "scater", "tcc"),
                         verbose = FALSE,
                         force = FALSE) {

  ## Generate paths to output folders
  datasetdir <- paste0(topdir, "/data-processed/", id)
  fastqcdir <- paste0(datasetdir, "/fastqc")
  salmondir <- paste0(datasetdir, "/salmon")
  umisdir <- paste0(datasetdir, "/umis")
  cutadaptdir <- paste0(datasetdir, "/cutadapt")
  scaterdir <- paste0(topdir, "/report-scater")
  maedir <- paste0(topdir, "/data-mae")
  multiqcdir <- paste0(topdir, "/report-multiqc")
  countsimqcdir <- paste0(topdir, "/report-countsimqc")
  tccdir <- paste0(topdir, "/data-tcc/", id)
  kallistodir <- paste0(tccdir, "/kallistotcc")
  
  ## Read run info downloaded from SRA
  x <- read.delim(paste0(topdir, "/data-raw/", id, "/", id, "_SraRunInfo.csv"), 
                  header = TRUE, as.is = TRUE, sep = ",")
  samples <- unique(x[, sncol])
  if(verbose) message( "Found ", length(samples), " samples to process.")


  any_updated <- 0
  for (smp in samples) {
    if(verbose) message( "Working on ", smp, " ..")
    ## Find all the SRA runs corresponding to this sample. They will be merged
    ## together in the analysis
    runs <- x$Run[x[, sncol] == smp]
    if(verbose) message( "Found ", length(runs), " runs.")
    
    ## Put together a file list
    if (rtype == "single") {
      files <- paste(paste0("<(./stream_ena ", runs, ".fastq)"), collapse = " ")
    } else if (rtype == "paired") {
      files1 <- paste(paste0("<(./stream_ena ", runs, "_1.fastq)"), collapse = " ")
      files2 <- paste(paste0("<(./stream_ena ", runs, "_2.fastq)"), collapse = " ")
      files <- list(f1 = files1, f2 = files2)
    } 
    if(verbose) message( "Files: ", files )
    
    if (any(c("fastqc", "salmon", "umis", "tcc") %in% aspects)) {
      if (force || 
          !(any(c(file.exists(paste0(fastqcdir, "/", smp, "/", smp, "_fastqc.html")),
                  file.exists(paste0(fastqcdir, "/", smp, "/", smp, "_1_fastqc.html")) && 
                  file.exists(paste0(fastqcdir, "/", smp, "/", smp, "_2_fastqc.html"))))) ||
          !(any(file.exists(c(paste0(salmondir, "/", smp, "/aux_info/meta_info.json"),
                              paste0(salmondir, "/", smp, "/aux/meta_info.json"),
                              paste0(umisdir, "/", smp, "/umi_counts.txt")))))) {
        any_updated <- 1
        
        ## Download fastq file(s) and save temporarily
        files <- download_fastq(rtype = rtype, outdir = tmpdir, smp = smp, files = files)
          
        ## Trim
        if (dotrim)
          files <- trim(rtype = rtype, cutadaptdir = cutadaptdir, smp = smp, 
                        adapterseq = adapterseq, cutadaptbin = cutadaptbin, 
                        fastqdir = tmpdir)
        
        ## FastQC
        if ("fastqc" %in% aspects)
          fastqc(rtype = rtype, fastqcdir = fastqcdir, smp = smp, 
                 files = files, fastqcbin = fastqcbin)
        
        ## Salmon
        if ("salmon" %in% aspects)
          quantify_salmon(rtype = rtype, files = files, 
                          salmondir = salmondir, smp = smp,
                          salmonbin = salmonbin, libtype = libtype, 
                          salmonindex = salmonindex, bias = bias)
        
        ## kallisto-tcc
        if ("tcc" %in% aspects) {
          if( verbose ) message("Running quantify_kallistotcc() on ", files)
          quantify_kallistotcc(rtype = rtype, files = files, 
                          kallistodir = kallistodir, smp = smp,
                          kallistobin = kallistobin, kallistoindex = kallistoindex)
        }
        
        ## RapMap + umis
        if ("umis" %in% aspects)
          quantify_umis(files = files, rapmapbin = rapmapbin, cell_barcodes = cell_barcodes, 
                        rapmapindex = rapmapindex, umis_transform = umis_transform, 
                        smp = smp, tmpdir = tmpdir, umisdir = umisdir)
      } else {
        message("Output files for ", smp, " already exist.")
      }
      sapply(files, function(f) 
        if (file.exists(f)) unlink(f)
      )
    }
  }
  
  if ("salmon" %in% aspects) {
    summary_table_salmon <- 
      compress_summarize_salmon(id = id, topdir = topdir, salmondir = salmondir,
                                datasetdir = datasetdir, any_updated = any_updated)
    write.table(cbind(sample = rownames(summary_table_salmon), summary_table_salmon), 
                file = paste0(datasetdir, "/summary_table_salmon.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }
  
  if ("umis" %in% aspects) {
    summary_table_rapmap <- 
      summarize_rapmap(id = id, umisdir = umisdir)
    write.table(summary_table_rapmap,
                file = paste0(datasetdir, "/summary_table_rapmap.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }

  ## ------------------------------------------------------------------------ ##
  ##                             MultiQC                                      ##
  ## ------------------------------------------------------------------------ ##
  if ("multiqc" %in% aspects) {
    message("Running MultiQC for ", id)
    mqc <- sprintf("bash -c '%s -o %s -n %s -f %s'",
                   multiqcbin, 
                   multiqcdir,
                   paste0(id, "_multiqc_report.html"),
                   datasetdir)
    system(mqc)
  }
  
  if ("mae" %in% aspects) {
    ## ------------------------------------------------------------------------ ##
    ##                  Generate MultiAssayExperiment                           ##
    ## ------------------------------------------------------------------------ ##
    if (dtype == "fl") {
      summary_table_salmon <- read.delim(paste0(datasetdir, "/summary_table_salmon.txt"),
                                         header = TRUE, as.is = TRUE)
      mae <- mae_tximport(id = id, salmondir = salmondir, topdir = topdir, 
                          txgenemap = txgenemap, phenofile = phenofile, geodata = geodata, 
                          gene_granges = gene_granges, tx_granges = tx_granges)
      mae@metadata <- list(genome = genome, 
                           organism = organism,
                           salmon_summary = summary_table_salmon,
                           creation_date = date())
    } else if (dtype == "umi") {
      summary_table_rapmap <- read.delim(paste0(datasetdir, "/summary_table_rapmap.txt"),
                                         header = TRUE, as.is = TRUE)
      mae <- mae_umis(id = id, umisdir = umisdir, topdir = topdir, 
                      txgenemap = txgenemap, phenofile = phenofile, geodata = geodata, 
                      gene_granges = gene_granges, tx_granges = tx_granges)
      mae@metadata <- list(genome = genome, 
                           organism = organism,
                           rapmap_summary = summary_table_rapmap, 
                           creation_date = date())
    }
    
    saveRDS(mae, file = paste0(maedir, "/", id, ".rds"))
    
    ## ------------------------------------------------------------------------ ##
    ##                  Write basic information to file                         ##
    ## ------------------------------------------------------------------------ ##
    infodf <- data.frame(nsamples = nrow(colData(mae)),
                         organism = organism,
                         genome = genome,
                         ntranscripts = nrow(assays(experiments(mae)[["tx"]])[["count"]]),
                         ngenes = nrow(assays(experiments(mae)[["gene"]])[["count"]]),
                         PMID = pmid,
                         datalink = datalink,
                         shortname = shortname,
                         description = description,
                         protocol = protocol,
                         protocoltype = protocoltype)
    write.table(t(infodf), file = paste0(datasetdir, "/dataset_info.txt"),
                row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
  
  ## ------------------------------------------------------------------------ ##
  ##                    Generate scater QC report                             ##
  ## ------------------------------------------------------------------------ ##
  if ("scater" %in% aspects) {
    mae <- readRDS(paste0(maedir, "/", id, ".rds"))
    message("Generating scater report for ", id)
    generate_report(id = id, maex = mae, phenoid = groupid, 
                    output_format = "html_document",
                    output_file = paste0(id, "_scater.html"),
                    output_dir = scaterdir,
                    nrw = nrw, lps = lps)
  }
  
  ## ------------------------------------------------------------------------ ##
  ##                    Generate countsimQC report                            ##
  ## ------------------------------------------------------------------------ ##
  if ("countsimQC" %in% aspects) {
    mae <- readRDS(paste0(maedir, "/", id, ".rds"))
    message("Generating countsimQC report for ", id)
    ddsList <- list(
      count_gene = DESeqDataSetFromMatrix(
        countData = round(assays(experiments(mae)[["gene"]])[["count"]]),
        colData = colData(mae),
        design = as.formula(paste0("~ ", paste(groupid, collapse = " + ")))),
      count_tx = DESeqDataSetFromMatrix(
        countData = round(assays(experiments(mae)[["tx"]])[["count"]]),
        colData = colData(mae),
        design = as.formula(paste0("~ ", paste(groupid, collapse = " + "))))
    )
    countsimQCReport(ddsList = ddsList, outputFile = paste0(id, "_countsimQC.html"),
                     outputDir = countsimqcdir, outputFormat = "html_document",
                     showCode = FALSE, forceOverwrite = TRUE, savePlot = FALSE,
                     description = id, subsampleSize = 250, maxNForDisp = 100)
  }
}

