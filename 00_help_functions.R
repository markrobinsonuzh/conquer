fastqc_single <- function(fastqc_dir, smp, files, fastqcbin, appd = "") {
  if (!file.exists(paste0(fastqc_dir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(fastqc_dir, "/", smp))
    system(mkd)
  }
  
  if (!file.exists(paste0(fastqc_dir, "/", smp, "/", smp, appd, "_fastqc.html"))) {
    fastqc <- sprintf("bash -c 'cat %s | %s --noextract -o %s -f fastq stdin:%s'",
                      files,
                      fastqcbin, 
                      paste0(fastqc_dir, "/", smp),
                      paste0(smp, appd))
    system(fastqc)
  } else {
    message("FastQC has already been run for ", smp, appd)
  }
}

salmon_single <- function(salmon_dir, smp, files, salmonbin, libtype, index) {
  if (!file.exists(paste0(salmon_dir, "/", smp, "/quant.sf"))) {
    salmon <- sprintf("bash -c '%s quant -p 10 -l %s -i %s -r <(cat %s) -o %s'",
                      salmonbin, 
                      libtype,
                      index,
                      files, 
                      paste0(salmon_dir, "/", smp))
    system(salmon)
  } else {
    message("Salmon has already been run for ", smp)
  }
}

salmon_paired <- function(salmon_dir, smp, files1, files2, salmonbin, libtype, index) {
  if (!file.exists(paste0(salmon_dir, "/", smp, "/quant.sf"))) {
    salmon <- sprintf("bash -c '%s quant -p 10 -l %s -i %s -1 <(cat %s) -2 <(cat %s) -o %s'",
                      salmonbin, 
                      libtype,
                      index,
                      files1,
                      files2,
                      paste0(salmon_dir, "/", smp))
    system(salmon)
  } else {
    message("Salmon has already been run for ", smp)
  }
}

fastqc_paired <- function(fastqc_dir, smp, files1, files2, fastqcbin) {
  fastqc_single(fastqc_dir = fastqc_dir, smp = smp, files = files1,
                fastqcbin = fastqcbin, appd = "_1")
  fastqc_single(fastqc_dir = fastqc_dir, smp = smp, files = files2,
                fastqcbin = fastqcbin, appd = "_2")
}

trim_single <- function(cutadapt_dir, smp, adapterseq, cutadaptbin, tmp_dir) {
  if (!file.exists(paste0(cutadapt_dir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(cutadapt_dir, "/", smp))
    system(mkd)
  }
  if (!file.exists(paste0(tmp_dir, "/", smp, ".trim.fastq"))) {
    ## Run trimming and save temporarily the resulting fastq file
    cutadapt <- sprintf("bash -c '%s -f fastq -m 15 -O 3 -a %s -o %s %s > %s'",
                        cutadaptbin, 
                        adapterseq,
                        paste0(tmp_dir, "/", smp, ".trim.fastq"),
                        paste0(tmp_dir, "/", smp, ".fastq"),
                        paste0(cutadapt_dir, "/", smp, "/", smp, "_cutadapt.txt"))
    system(cutadapt)
    
    ## Remove the original fastq file
    unlink(paste0(tmp_dir, "/", smp, ".fastq"))
  }
}

trim_paired <- function(cutadapt_dir, smp, adapterseq, cutadaptbin, tmp_dir) {
  if (!file.exists(paste0(cutadapt_dir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(cutadapt_dir, "/", smp))
    system(mkd)
  }
  
  if (!all(file.exists(paste0(tmp_dir, "/", smp, c("_1", "_2"), ".fastq")))) {
    ## Run trimming adn same temporarily the resulting fastq files
    cutadapt <- sprintf("bash -c '%s -f fastq -m 15 -O 3 -a %s -A %s -o %s -p %s %s %s > %s'",
                        cutadaptbin, 
                        adapterseq,
                        adapterseq,
                        paste0(tmp_dir, "/", smp, "_1.trim.fastq"),
                        paste0(tmp_dir, "/", smp, "_2.trim.fastq"),
                        paste0(tmp_dir, "/", smp, "_1.fastq"),
                        paste0(tmp_dir, "/", smp, "_2.fastq"),
                        paste0(cutadapt_dir, "/", smp, "/", smp, "_cutadapt.txt"))
    system(cutadapt)
    
    ## Remove the original fastq files
    unlink(paste0(tmp_dir, "/", smp, "_1.fastq"))
    unlink(paste0(tmp_dir, "/", smp, "_2.fastq"))
  }
}

