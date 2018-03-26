## -------------------------------------------------------------------------- ##
##                            Download fastq(s)                               ##
## -------------------------------------------------------------------------- ##
#' Download fastq file(s) and save in temporary directory
#' 
#' The name of the downloaded file(s) will be smp.fastq, where smp is the
#' provided Sample ID. If rtype = "paired", _1 and _2 will be appended to the
#' sample ID.
#' 
#' @param rtype "single" or "paired"
#' @param outdir The directory where the files are saved
#' @param smp Sample ID
#' @param files String (or vector of strings) that point to the downloaded
#'   files. Typically, of the form "<(./stream_ena GSM12345.fasta)"
#'
#' @examples
#' ## Download the file ERR1042832.fastq from SRA and save it as
#' ## tmp/E7.8.350.fastq
#' download_fastq(rtype = "single", outdir = "tmp", smp = "E7.8.350", 
#'                files = "<(./stream_ena ERR1042832.fastq)") 
#'
#' @return The name(s) of the downloaded fastq file(s).
download_fastq <- function(rtype, outdir, smp, files) {
  if (rtype == "single") {
    message("Downloading fastq file for ", smp)
    dwl <- sprintf("bash -c 'cat %s > %s'",
                   files,
                   paste0(outdir, "/", smp, ".fastq"))
    system(dwl)
    return(paste0(outdir, "/", smp, ".fastq"))
  } else if (rtype == "paired") {
    message("Downloading fastq files for ", smp)
    dwl1 <- sprintf("bash -c 'cat %s > %s'",
                    files$f1,
                    paste0(outdir, "/", smp, "_1.fastq"))
    system(dwl1)
    dwl2 <- sprintf("bash -c 'cat %s > %s'",
                    files$f2,
                    paste0(outdir, "/", smp, "_2.fastq"))
    system(dwl2)
    return(list(f1 = paste0(outdir, "/", smp, "_1.fastq"), 
                f2 = paste0(outdir, "/", smp, "_2.fastq")))
  }
}

## -------------------------------------------------------------------------- ##
##                                Run kallisto TCC                            ##
## -------------------------------------------------------------------------- ##
#' Compute equivalence classes and quantify abundances w/ kallisto
#' 
#' @param rtype "single" or "paired"
#' @param files Names of the fastq files to use for the quantification
#' @param smp Sample ID
#' @param kallistobin Path to kallisto binary
#' @param kallistoindex Path to kallisto index
#' 
#' @return Returns nothing, but generates kallisto output files in the
#'   kallistodir/smp directory.
#'   
quantify_kallistotcc <- function(rtype, files, kallistodir, smp, kallistobin,
                            kallistoindex) {
  kallisto_files <- c("pseudoalignments.ec", "pseudoalignments.tsv", "run_info.json")
  if ( !file.exists(c(paste0(kallistodir, "/", smp, "/run_info.json"))) &
       !file.exists(c(paste0(kallistodir, "/", smp, "/run_info.json.gz"))) ) {
    out_dir <- paste0(kallistodir, "/", smp)
    if (rtype == "single") {
      kallisto <- sprintf("bash -c '%s pseudo -t 8 --single -l 200 -s 30 -i %s -o %s %s'",
                        kallistobin, kallistoindex, out_dir, files)
    } else if (rtype == "paired") {
      kallisto <- sprintf("bash -c '%s pseudo -t 10 -i %s -o %s %s %s'",
                        kallistobin, kallistoindex, out_dir,
                        files$f1, files$f2)
    }
    cat(kallisto)
    system(kallisto)
    # compress files
    kallisto_compress <- paste0("gzip ", paste(file.path(out_dir,kallisto_files), collapse=" "))
    system(kallisto_compress)
  } else {
    message("kallisto has already been run for ", smp)
  }
}


## -------------------------------------------------------------------------- ##
##                                Compile kallisto TCC                        ##
## -------------------------------------------------------------------------- ##
#' Compile a set of kallisto .ec and .tsv and 
#'
#' @param base_dir directory where subdirectories of data (that contain .ec.gz, .tsv.gz) are
#' @param mae_dir directory with the existing salmon-based MultiAssayExperiment objects are (data-mae)
#' @param out_dir directory to deposit the TCC SummarizedExperiment object to
#' @param experiment_id identifier of the experiment
#' @param verbose boolean, whether to write message() statements as the program progresses
#'
#' @return returns invisibly the SummarizeExperiment object, generates RDS file of the 
#'   object into the out_dir directory
#'

compile_tcc_counts <- function(base_dir, mae_dir, out_dir, experiment_id, verbose=TRUE) {


  # get filenames
  ec_files <- dir(base_dir, "pseudoalignments.ec.gz", recursive=TRUE, full.names=TRUE)
  count_files <- dir(base_dir, "pseudoalignments.tsv.gz", recursive=TRUE, full.names=TRUE)

  stopifnot(length(ec_files)==length(count_files))

  if(verbose) message( "Found ", length(ec_files), " samples.")

  # note this 8 below is because the sample id is after the 11th "/"
  # /home/Shared/data/seq/conquer/database/data-tcc/SRP073808/kallistotcc/SRR3952971/..
  extract_delim <- function(x, delim="/", ind=11) {
    sapply( strsplit(x,delim), .subset, ind)
  }

  sample_id_ec <- extract_delim(ec_files, "/", 11)
  sample_id_cnt <- extract_delim(count_files, "/", 11)

  stopifnot(sample_id_ec==sample_id_cnt)

  ec_col <- cols(
    label = col_integer(),
    trans_in_class = col_character()
  )

  counts_col <- cols(
    label = col_integer(),
    count = col_integer()
  )


  if(verbose) message( "Reading ", length(ec_files), " samples to process.")

  dfs <- mapply( function(u,v,z) {
    cat(".")
    ec <- readr::read_tsv(u, col_names=c("label", "trans_in_class"), 
                          col_types = ec_col, progress=FALSE )
    cnt <- readr::read_tsv(v, col_names=c("label", "count"),
                           col_types = counts_col, progress=FALSE )
    k <- cnt$count > 0
    df <- data.frame(ec=ec$trans_in_class[k], 
                     counts=cnt$count[k], stringsAsFactors=FALSE)
    colnames(df)[2] <- z
    df
  }, ec_files, count_files, sample_id_ec, SIMPLIFY=FALSE)

  if(verbose) message( "Performing full_join().")
  df_merged <- dfs %>%
                 Reduce(function(df1,df2) full_join(df1,df2,by="ec"), .) %>%
                 column_to_rownames("ec") %>%
                 replace(., is.na(.), 0)

  mae_rds <- paste0(experiment_id, ".rds")
  if(verbose) message( "Reading ", file.path(mae_dir, mae_rds), ".")
  mae <- readRDS( file.path(mae_dir, mae_rds) )

  if(verbose) message( "Constructing SummarizedExperiment.")
  samples <- intersect(colnames(df_merged), mae$Run)

  md <- metadata(mae)

  m1 <- match(samples , colnames(df_merged))
  m2 <- match(samples, mae$Run)

  cd <- colData(mae)
  x <- as.matrix(df_merged)

  se <- SummarizedExperiment(assays=SimpleList(tcc=x[,m1]),
                             colData=cd[m2,])

  tcc_rds <- paste0(experiment_id, ".rds")
  if(verbose) message( "Writing ", file.path(out_dir, tcc_rds), ".")
  saveRDS(se, file.path(out_dir, tcc_rds))
  invisible(se)
}




## -------------------------------------------------------------------------- ##
##                              Trim adapters                                 ##
## -------------------------------------------------------------------------- ##
#' Trim adapter sequences using cutadapt
#' 
#' @param rtype "single" or "paired"
#' @param cutadaptdir Directory to write trimming information to
#' @param smp Sample ID
#' @param adapterseq Adapter sequence
#' @param cutadaptbin Path to the cutadapt binary
#' @param fastqdir Directory where fastq files to be trimmed are located. The
#'   name(s) of these fastq files should be of the form smp.fastq or
#'   smp_1.fastq/smp_2.fastq.
#'
#' @return The names of the trimmed fastq file(s). Also write trimming
#'   information to the cutadaptdir directory.
#' 
trim <- function(rtype, cutadaptdir, smp, adapterseq, cutadaptbin, fastqdir) {
  message("Trimming fastq file(s) for ", smp)
  if (!file.exists(paste0(cutadaptdir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(cutadaptdir, "/", smp))
    system(mkd)
  }
  if (rtype == "single") {
    if (!file.exists(paste0(fastqdir, "/", smp, ".trim.fastq"))) {
      ## Run trimming and save temporarily the resulting fastq file
      cutadapt <- sprintf("bash -c '%s -f fastq -m 15 -O 3 -a %s -o %s %s > %s'",
                          cutadaptbin, 
                          adapterseq,
                          paste0(fastqdir, "/", smp, ".trim.fastq"),
                          paste0(fastqdir, "/", smp, ".fastq"),
                          paste0(cutadaptdir, "/", smp, "/", smp, "_cutadapt.txt"))
      system(cutadapt)
      ## Remove the original fastq file
      unlink(paste0(fastqdir, "/", smp, ".fastq"))
    }
    files <- paste0(fastqdir, "/", smp, ".trim.fastq")
  } else if (rtype == "paired") {
    if (!all(file.exists(paste0(fastqdir, "/", smp, c("_1", "_2"), ".trim.fastq")))) {
      ## Run trimming and save temporarily the resulting fastq files
      cutadapt <- sprintf("bash -c '%s -f fastq -m 15 -O 3 -a %s -A %s -o %s -p %s %s %s > %s'",
                          cutadaptbin, 
                          adapterseq,
                          adapterseq,
                          paste0(fastqdir, "/", smp, "_1.trim.fastq"),
                          paste0(fastqdir, "/", smp, "_2.trim.fastq"),
                          paste0(fastqdir, "/", smp, "_1.fastq"),
                          paste0(fastqdir, "/", smp, "_2.fastq"),
                          paste0(cutadaptdir, "/", smp, "/", smp, "_cutadapt.txt"))
      system(cutadapt)
      ## Remove the original fastq files
      unlink(paste0(fastqdir, "/", smp, "_1.fastq"))
      unlink(paste0(fastqdir, "/", smp, "_2.fastq"))
    }
    files = list(f1 = paste0(fastqdir, "/", smp, "_1.trim.fastq"),
                 f2 = paste0(fastqdir, "/", smp, "_2.trim.fastq"))
  }
  files
}

## -------------------------------------------------------------------------- ##
##                                Run FastQC                                  ##
## -------------------------------------------------------------------------- ##
fastqc_single <- function(fastqcdir, smp, files, fastqcbin, appd = "") {
  if (!file.exists(paste0(fastqcdir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(fastqcdir, "/", smp))
    system(mkd)
  }
  
  if (!file.exists(paste0(fastqcdir, "/", smp, "/", smp, appd, "_fastqc.html"))) {
    fastqc <- sprintf("bash -c 'cat %s | %s --noextract -o %s -t 10 -f fastq stdin:%s'",
                      files,
                      fastqcbin, 
                      paste0(fastqcdir, "/", smp),
                      paste0(smp, appd))
    system(fastqc)
  } else {
    message("FastQC has already been run for ", smp, appd)
  }
}

#' Run FastQC on fastq file(s)
#' 
#' @param rtype "single" or "paired"
#' @param fastqcdir Directory where FastQC reports should be written
#' @param smp Sample ID
#' @param files Name(s) of fastq file(s) to run FastQC on
#' @param fastqcbin Path to FastQC binary
#' @param appd String to append to report name
#' 
#' @param return Nothing is returned, but FastQC reports are generated in the
#'   fastqc_dir directory.
#'   
fastqc <- function(rtype, fastqcdir, smp, files, fastqcbin) {
  if (rtype == "single") {
    fastqc_single(fastqcdir, smp, files, fastqcbin, appd = "")
  } else if (rtype == "paired") {
    fastqc_single(fastqcdir = fastqcdir, smp = smp, files = files$f1,
                  fastqcbin = fastqcbin, appd = "_1")
    fastqc_single(fastqcdir = fastqcdir, smp = smp, files = files$f2,
                  fastqcbin = fastqcbin, appd = "_2")
  }
}

## -------------------------------------------------------------------------- ##
##                                Run Salmon                                  ##
## -------------------------------------------------------------------------- ##
#' Quantify transcript abundance with Salmon
#' 
#' @param rtype "single" or "paired"
#' @param files Names of the fastq files to use for the quantification
#' @param smp Sample ID
#' @param salmonbin Path to Salmon binary
#' @param libtype The \code{LIBTYPE} argument passed to Salmon
#' @param index Path to Salmon index
#' @param bias Whether or not to use the --seqBias argument of Salmon  
#' 
#' @return Returns nothing, but generates Salmon output files in the
#'   salmondir/smp directory.
#'   
quantify_salmon <- function(rtype, files, salmondir, smp, salmonbin, 
                            libtype, salmonindex, bias = FALSE) {
  if (!any(file.exists(c(paste0(salmondir, "/", smp, "/aux_info/meta_info.json"),
                         paste0(salmondir, "/", smp, "/aux/meta_info.json"))))) {
    if (rtype == "single") {
      salmon <- sprintf("bash -c '%s quant -p 10 -l %s -i %s -r <(cat %s) -o %s %s'",
                        salmonbin, 
                        libtype,
                        salmonindex,
                        files, 
                        paste0(salmondir, "/", smp),
                        ifelse(bias, "--seqBias", ""))
      system(salmon)
    } else if (rtype == "paired") {
      salmon <- sprintf("bash -c '%s quant -p 10 -l %s -i %s -1 <(cat %s) -2 <(cat %s) -o %s %s'",
                        salmonbin, 
                        libtype,
                        salmonindex,
                        files$f1,
                        files$f2,
                        paste0(salmondir, "/", smp),
                        ifelse(bias, "--seqBias", ""))
      system(salmon)
    } 
  } else {
    message("Salmon has already been run for ", smp)
  }
}

## -------------------------------------------------------------------------- ##
##                 Compress and summarize Salmon output                       ##
## -------------------------------------------------------------------------- ##
compress_summarize_salmon <- function(id, topdir, salmondir, datasetdir,
                                      any_updated = 1) {
  ## Compress all Salmon output in a tar archive
  if (any_updated == 1) {
    message("Compressing Salmon output for ", id)
    targz <- sprintf("bash -c 'tar -C %s/data-processed/ -czf %s %s'",
                     topdir, 
                     paste0(datasetdir, "/", id, "_salmon.tar.gz"),
                     paste0(id, "/salmon"))
    system(targz)
  }
  
  ## Create summary table from Salmon parameters and results
  message("Creating Salmon summary table for ", id)
  smps <- list.files(salmondir, full.names = TRUE)
  names(smps) <- basename(smps)
  summary_table_salmon <- as.data.frame(t(sapply(smps, function(s) {
    cmdinfo <- fromJSON(file = paste0(s, "/cmd_info.json"))
    cmdinfo <- c(salmon_version = cmdinfo[["salmon_version"]],
                 libtype = cmdinfo[["libType"]],
                 index = basename(cmdinfo[["index"]]),
                 seqBias = ifelse("seqBias" %in% names(cmdinfo), "TRUE", "FALSE"))
    if (file.exists(paste0(s, "/aux_info"))) {
      metainfo <- fromJSON(file = paste0(s, "/aux_info/meta_info.json"))
    } else {
      metainfo <- fromJSON(file = paste0(s, "/aux/meta_info.json"))
    }
    metainfo <- c(num_processed = format(metainfo[["num_processed"]], scientific = FALSE),
                  num_mapped = format(metainfo[["num_mapped"]], scientific = FALSE),
                  percent_mapped = round(metainfo[["percent_mapped"]], 3))
    c(cmdinfo, metainfo)
  })))
  summary_table_salmon
}

## -------------------------------------------------------------------------- ##
##         Generate MultiAssayExperiment object from Salmon results           ##
## -------------------------------------------------------------------------- ##
## Make sure that column type is determined from all values in a column when
## importing Salmon files
read_tsv2 <- function(...) readr::read_tsv(..., guess_max = 100000, progress = FALSE,
                                           col_types = list(
                                             Name = col_character(),
                                             Length = col_integer(),
                                             EffectiveLength = col_double(),
                                             TPM = col_double(),
                                             NumReads = col_double()
                                           ))

mae_tximport <- function(id, salmondir, topdir, txgenemap, geodata,
                         phenofile, gene_granges, tx_granges) {
  message("Reading expression levels for ", id)
  files <- paste0(list.files(salmondir, full.names = TRUE), "/quant.sf")
  names(files) <- basename(gsub("/quant.sf", "", files))
  txi_tx <- tximport(files = files, type = "salmon", txIn = TRUE, txOut = TRUE, 
                     dropInfReps = TRUE, importer = read_tsv2)
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
                              colData = droplevels(meta))
  mae
}

## -------------------------------------------------------------------------- ##
##          Generate MultiAssayExperiment object from umis results            ##
## -------------------------------------------------------------------------- ##
mae_umis <- function(id, umisdir, topdir, txgenemap, geodata,
                     phenofile, gene_granges, tx_granges) {
  message("Reading expression levels for ", id)
  files <- paste0(list.files(umisdir, full.names = TRUE), "/umi_counts.txt")
  names(files) <- basename(gsub("/umi_counts.txt", "", files))
  
  tmp <- lapply(files, function(f) {
    read.delim(f, sep = ",", as.is = TRUE, header = TRUE)
  })
  for (nm in names(tmp)) {
    w <- grep(nm, colnames(tmp[[nm]]), invert = TRUE)
    colnames(tmp[[nm]])[w] <- paste0(nm, "_", colnames(tmp[[nm]])[w])
    colnames(tmp[[nm]])[grep("gene", colnames(tmp[[nm]]))] <- "tx"
  }
  txi_tx <- Reduce(function(...) dplyr::full_join(..., by = "tx"), tmp
  ) %>% dplyr::left_join(txgenemap)
  txi_gene <- txi_tx %>% dplyr::select(-tx) %>% group_by(gene) %>% 
    dplyr::summarise_all(funs(sum)) %>% as.data.frame()
  txi_tx <- txi_tx %>% dplyr::select(-gene)
  rownames(txi_tx) <- txi_tx$tx
  txi_tx$tx <- NULL
  rownames(txi_gene) <- txi_gene$gene
  txi_gene$gene <- NULL

  if (geodata == TRUE) {
    geo <- getGEO(filename = paste0(topdir, "/data-raw/", id, "/", id, "_series_matrix.txt.gz"),
                  getGPL = FALSE)
    meta <- pData(geo)
  } else {
    meta <- read.delim(phenofile, header = TRUE, row.names = 1, as.is = TRUE)
  }
  stopifnot(all(sapply(strsplit(colnames(txi_tx), "_"), .subset, 1) %in% rownames(meta)))
  stopifnot(all(colnames(txi_tx) == colnames(txi_gene)))
  meta <- meta[match(sapply(strsplit(colnames(txi_tx), "_"), .subset, 1), rownames(meta)), ]
  rownames(meta) <- colnames(txi_tx)
  
  generse <- SummarizedExperiment(assays = list(count = as.matrix(txi_gene)),
                                  rowRanges = gene_granges[rownames(txi_gene)])
  
  txrse <- SummarizedExperiment(assays = list(count = as.matrix(txi_tx)),
                                rowRanges = tx_granges[rownames(txi_tx)])
  
  ## Generate MultiAssayExperiment
  mae <- MultiAssayExperiment(experiments = list(gene = generse,
                                                 tx = txrse),
                              colData = droplevels(meta))
  mae
}

## -------------------------------------------------------------------------- ##
##                        Generate scater report                              ##
## -------------------------------------------------------------------------- ##
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
