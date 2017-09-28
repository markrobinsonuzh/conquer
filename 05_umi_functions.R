quantify_umis <- function(files, umisdir, rapmapbin, rapmapindex, 
                          umis_transform, smp, tmpdir, cell_barcodes) {
  ## umis fastqtransform
  if (length(grep("CEL-Seq2/transform.json", umis_transform)) != 0) {
    stopifnot(length(files) == 1)
    cmd <- sprintf("umis fastqtransform --cores 10 --demuxed_cb %s %s %s > %s",
                   gsub("\\.fastq$|\\.fq$", "", basename(files)),
                        umis_transform,
                        files,
                        paste0(tmpdir, "/", smp, "_transformed.fastq")
    )
    message(cmd)
    system(cmd)
  } else if (length(grep("CEL-Seq/transform.json", umis_transform)) != 0) {
    stopifnot(length(files) == 2)
    cmd <- sprintf("bash -c 'umis fastqtransform --cores 10 %s %s %s > %s'",
                   umis_transform,
                   files$f1,
                   files$f2,
                   paste0(tmpdir, "/", smp, "_transformed1.fastq")
    )
    message(cmd)
    system(cmd)
    cmd <- sprintf("umis cb_filter --bc1 %s --nedit 1 --cores 10 %s > %s",
                   cell_barcodes, 
                   paste0(tmpdir, "/", smp, "_transformed1.fastq"),
                   paste0(tmpdir, "/", smp, "_transformed.fastq")
    )
    message(cmd)
    system(cmd)
    unlink(paste0(tmpdir, "/", smp, "_transformed1.fastq"))
  }
  
  
  cmd <- sprintf("%s quasimap -t 10 -o %s -r %s -i %s",
                 rapmapbin,
                 paste0(tmpdir, "/", smp, "_rapmap.sam"),
                 paste0(tmpdir, "/", smp, "_transformed.fastq"),
                 rapmapindex
  )
  message(cmd)
  system(cmd)
  
  ## Read sam file header and extract RapMap version number
  cmd <- sprintf("samtools view -S -H %s > %s",
                 paste0(tmpdir, "/", smp, "_rapmap.sam"),
                 paste0(tmpdir, "/", smp, "_rapmap_header.txt"))
  message(cmd)
  system(cmd)
  hdr <- readLines(paste0(tmpdir, "/", smp, "_rapmap_header.txt"))
  tmp <- strsplit(hdr[grep("@PG", hdr)], "\t")[[1]]
  rapmap_version <- gsub("VN:", "", tmp[grep("VN:", tmp)])
  
  cmd <- sprintf("mkdir -p %s",
                 paste0(umisdir, "/", smp)
  )
  system(cmd)
  
  cmd <- sprintf("umis tagcount %s %s",
                 paste0(tmpdir, "/", smp, "_rapmap.sam"),
                 paste0(umisdir, "/", smp, "/umi_counts.txt")
  )
  message(cmd)
  system(cmd)
  
  unlink(paste0(tmpdir, "/", smp, "_rapmap.sam"))
  unlink(paste0(tmpdir, "/", smp, "_rapmap_header.txt"))
  unlink(paste0(tmpdir, "/", smp, "_transformed.fastq"))
  
  if (is.null(cell_barcodes)) cell_barcodes <- "None"
  
  write.table(data.frame(sample = smp, rapmap_version = rapmap_version, 
                         rapmap_index = rapmapindex, 
                         umis_transform = umis_transform, 
                         cell_barcodes = basename(cell_barcodes)),
              file = paste0(umisdir, "/", smp, "/rapmap_settings.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
}

summarize_rapmap <- function(id = id, umisdir = umisdir) {
  message("Creating RapMap summary table for ", id)
  smps <- list.files(umisdir, full.names = TRUE)
  names(smps) <- basename(smps)
  summary_table_rapmap <- data.frame(do.call(rbind, lapply(smps, function(s) {
    read.delim(paste0(s, "/rapmap_settings.txt"), header = TRUE, as.is = TRUE)
  })), stringsAsFactors = FALSE)
  summary_table_rapmap$rapmap_index <- basename(summary_table_rapmap$rapmap_index)
  summary_table_rapmap
  
}