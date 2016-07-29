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
    salmon <- sprintf("bash -c '%s quant -l %s -i %s -r <(cat %s) -o %s'",
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
    salmon <- sprintf("bash -c '%s quant -l %s -i %s -1 <(cat %s) -2 <(cat %s) -o %s'",
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
    # ## Download fastq file and save temporarily
    # dwl <- sprintf("bash -c 'cat %s > %s'",
    #                files,
    #                paste0(tmp_dir, "/", smp, ".fastq"))
    # system(dwl)
    
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
    
    # ## Modify the filename in the cutadapt output file
    # chf <- sprintf("bash -c 'sed -i \'s:/dev/fd/63:%s:\' %s'",
    #                paste0(smp, ".fastq"),
    #                paste0(cutadapt_dir, "/", smp, "/", smp, "_cutadapt.txt"))
    # system(chf)
  }
}

trim_paired <- function(cutadapt_dir, smp, adapterseq, cutadaptbin, tmp_dir) {
  if (!file.exists(paste0(cutadapt_dir, "/", smp))) {
    mkd <- sprintf("mkdir -p %s", paste0(cutadapt_dir, "/", smp))
    system(mkd)
  }
  
  if (!all(file.exists(paste0(tmp_dir, "/", smp, c("_1", "_2"), ".fastq")))) {
    # ## Download fastq files and save temporarily
    # dwl1 <- sprintf("bash -c 'cat %s > %s'",
    #                 files1,
    #                 paste0(tmp_dir, "/", smp, "_1.fastq"))
    # system(dwl1)
    # dwl2 <- sprintf("bash -c 'cat %s > %s'",
    #                 files2,
    #                 paste0(tmp_dir, "/", smp, "_2.fastq"))
    # system(dwl2)
    
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
    
    # ## Modify the filename in the cutadapt output file
    # chf <- sprintf("bash -c 'sed -i \'s:/dev/fd/63:%s:\' %s'",
    #                paste0(smp, "_1.fastq"),
    #                paste0(cutadapt_dir, "/", smp, "/", smp, "_cutadapt.txt"))
    # system(chf)
    # ## Modify the filename in the cutadapt output file
    # chf <- sprintf("bash -c 'sed -i \'s:/dev/fd/62:%s:\' %s'",
    #                paste0(smp, "_2.fastq"),
    #                paste0(cutadapt_dir, "/", smp, "/", smp, "_cutadapt.txt"))
    # system(chf)
  }
}

