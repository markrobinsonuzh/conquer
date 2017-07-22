## The code below builds the Salmon and RapMap indexes that will be used for 
## mapping and quantification. It also generates transcript-to-gene maps as well
## as GRanges objects characterizing the genes and transcripts

## Define reference files and corresponding Ensembl version
human_cdna_fa <- "reference-files/homo-sapiens/Homo_sapiens.GRCh38.cdna.all.fa"
human_ncrna_fa <- "reference-files/homo-sapiens/Homo_sapiens.GRCh38.ncrna.fa"
mouse_cdna_fa <- "reference-files/mus-musculus/Mus_musculus.GRCm38.cdna.all.fa"
mouse_ncrna_fa <- "reference-files/mus-musculus/Mus_musculus.GRCm38.ncrna.fa"
zebrafish_cdna_fa <- "reference-files/danio-rerio/Danio_rerio.GRCz10.cdna.all.fa"
zebrafish_ncrna_fa <- "reference-files/danio-rerio/Danio_rerio.GRCz10.ncrna.fa"
ercc_fa <- "reference-files/ERCC92/ERCC92.fa"
ercc_gtf <- "reference-files/ERCC92/ERCC92.gtf"
human_ensembl_version <- 84
mouse_ensembl_version <- 84
zebrafish_ensembl_version <- 87

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(rtracklayer))

## -------------------------------------------------------------------------- ##
##                            Build RapMap index                              ##
## -------------------------------------------------------------------------- ##

## Path to the RapMap binary
rapmapbin <- "/home/charlotte/software/RapMap/bin/rapmap"
rapmapversion <- "0.5.0"

## Human cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 31 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", human_cdna_fa, " ", human_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(human_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", rapmapversion, 
                                           ".ridx"), human_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 15 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", human_cdna_fa, " ", human_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(human_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", rapmapversion, 
                                           ".ridx"), human_cdna_fa))
system(cmd)

## Mouse cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 31 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", mouse_cdna_fa, " ", mouse_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", rapmapversion, 
                                           ".ridx"), mouse_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 15 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", mouse_cdna_fa, " ", mouse_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", rapmapversion, 
                                           ".ridx"), mouse_cdna_fa))
system(cmd)

## Zebrafish cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 31 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", zebrafish_cdna_fa, " ", zebrafish_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", rapmapversion, 
                                           ".ridx"), zebrafish_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s quasiindex -x 5 -k 15 -t %s -i %s'",
               rapmapbin,
               paste0("<(cat ", zebrafish_cdna_fa, " ", zebrafish_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", rapmapversion, 
                                           ".ridx"), zebrafish_cdna_fa))
system(cmd)

## -------------------------------------------------------------------------- ##
##                            Build Salmon index                              ##
## -------------------------------------------------------------------------- ##

## Path to the Salmon binary
#salmonbin <- "/usr/local/software/SalmonBeta-0.6.1_DebianSqueeze/bin/salmon"
#salmonbin <- "software/Salmon-0.7.2_linux_x86_64/bin/salmon"
#salmonversion <- "0.7.2"
salmonbin <- "/home/charlotte/software/Salmon-0.8.2_linux_x86_64/bin/salmon"
salmonversion <- "0.8.2"

## Human cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi'",
               salmonbin,
               paste0("<(cat ", human_cdna_fa, " ", human_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(human_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", salmonversion, 
                                           ".sidx"), human_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi -k 15'",
               salmonbin,
               paste0("<(cat ", human_cdna_fa, " ", human_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(human_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", salmonversion, 
                                           ".sidx"), human_cdna_fa))
system(cmd)


## Mouse cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi'",
               salmonbin,
               paste0("<(cat ", mouse_cdna_fa, " ", mouse_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", salmonversion,
                                           ".sidx"), mouse_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi -k 15'",
               salmonbin,
               paste0("<(cat ", mouse_cdna_fa, " ", mouse_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", salmonversion, 
                                           ".sidx"), mouse_cdna_fa))
system(cmd)

## Zebrafish cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi'",
               salmonbin,
               paste0("<(cat ", zebrafish_cdna_fa, " ", zebrafish_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.", salmonversion, 
                                           ".sidx"), zebrafish_cdna_fa))
system(cmd)

## Index for shorter reads
cmd <- sprintf("bash -c '%s index -t %s -i %s --type quasi -k 15'",
               salmonbin,
               paste0("<(cat ", zebrafish_cdna_fa, " ", zebrafish_ncrna_fa, " ", ercc_fa, ")"),
               gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.k15.", salmonversion,
                                           ".sidx"), zebrafish_cdna_fa))
system(cmd)

## -------------------------------------------------------------------------- ##
##                          Generate tx-to-gene map                           ##
## -------------------------------------------------------------------------- ##

ercc <- readDNAStringSet(ercc_fa)
ercc <- data.frame(tx = names(ercc), gene = names(ercc), stringsAsFactors = FALSE)

## Human
cdna <- readDNAStringSet(human_cdna_fa)
ncrna <- readDNAStringSet(human_ncrna_fa)
cdna <- data.frame(t(sapply(as.character(names(cdna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
ncrna <- data.frame(t(sapply(as.character(names(ncrna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
txgenemap <- rbind(cdna, ncrna, ercc)
saveRDS(txgenemap, 
        file = gsub("cdna.all.fa$", paste0(human_ensembl_version,
                                           ".cdna.ncrna.ercc92.txgenemap.rds"), human_cdna_fa))

## Mouse
cdna <- readDNAStringSet(mouse_cdna_fa)
ncrna <- readDNAStringSet(mouse_ncrna_fa)
cdna <- data.frame(t(sapply(as.character(names(cdna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
ncrna <- data.frame(t(sapply(as.character(names(ncrna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
txgenemap <- rbind(cdna, ncrna, ercc)
saveRDS(txgenemap, 
        file = gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.txgenemap.rds"), mouse_cdna_fa))

## Zebrafish
cdna <- readDNAStringSet(zebrafish_cdna_fa)
ncrna <- readDNAStringSet(zebrafish_ncrna_fa)
cdna <- data.frame(t(sapply(as.character(names(cdna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
ncrna <- data.frame(t(sapply(as.character(names(ncrna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), gene = ifelse(length(gene) != 0, gene, NA))
})), stringsAsFactors = FALSE)
txgenemap <- rbind(cdna, ncrna, ercc)
saveRDS(txgenemap, 
        file = gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.txgenemap.rds"), zebrafish_cdna_fa))

## -------------------------------------------------------------------------- ##
##                          Generate GRanges object                           ##
## -------------------------------------------------------------------------- ##
## Human
cdna <- readDNAStringSet(human_cdna_fa)
ncrna <- readDNAStringSet(human_ncrna_fa)

nm <- c(names(cdna), names(ncrna))
info <- data.frame(t(sapply(nm, function(w) {
  w <- strsplit(w, " ")[[1]]
  transcript <- w[1]
  gene <- gsub("gene:", "", w[grep("^gene:", w)])
  position <- gsub("chromosome:", "", w[grep("^chromosome:", w)])
  if (length(position) == 0) position <- gsub("scaffold:", "", w[grep("^scaffold:", w)])
  symbol <- gsub("gene_symbol:", "", w[grep("^gene_symbol", w)])
  c(transcript = transcript, 
    gene = gene,
    genome = strsplit(position, ":")[[1]][1], 
    chromosome = strsplit(position, ":")[[1]][2],
    start = strsplit(position, ":")[[1]][3],
    end = strsplit(position, ":")[[1]][4],
    strand = ifelse(strsplit(position, ":")[[1]][5] == "1", "+", "-"),
    symbol = symbol
  )
})), stringsAsFactors = FALSE)
rownames(info) <- NULL
info$start <- as.numeric(info$start)
info$end <- as.numeric(info$end)

txgr <- GRanges(seqnames = info$chromosome, 
                ranges = IRanges(start = info$start, end = info$end), 
                strand = info$strand)
mcols(txgr) <- info[, c("transcript", "gene", "genome", "symbol")]

erccgtf <- import(ercc_gtf, format = "gtf")
ercctx <- erccgtf
mcols(ercctx) <- data.frame(transcript = ercctx$gene_id,
                            gene = ercctx$gene_id,
                            genome = ercctx$source,
                            symbol = ercctx$gene_id)

txgr <- suppressWarnings(c(txgr, ercctx))
names(txgr) <- txgr$transcript

geneinfo <- info %>% group_by(gene) %>% 
  summarise(genome = unique(genome),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            strand = unique(strand),
            symbol = unique(symbol)) %>% 
  ungroup()

ggr <- GRanges(seqnames = geneinfo$chromosome,
               ranges = IRanges(start = geneinfo$start, end = geneinfo$end),
               strand = geneinfo$strand)
mcols(ggr) <- geneinfo[, c("gene", "genome", "symbol")]

erccgene <- erccgtf
mcols(erccgene) <- data.frame(gene = erccgene$gene_id,
                              genome = erccgene$source,
                              symbol = erccgene$gene_id)
ggr <- suppressWarnings(c(ggr, erccgene))
names(ggr) <- ggr$gene

gene_granges <- ggr
tx_granges <- txgr
saveRDS(list(gene_granges = gene_granges, tx_granges = tx_granges),
        file = gsub("cdna.all.fa$", paste0(human_ensembl_version, 
                                           ".cdna.ncrna.ercc92.granges.rds"), human_cdna_fa))

## Mouse
cdna <- readDNAStringSet(mouse_cdna_fa)
ncrna <- readDNAStringSet(mouse_ncrna_fa)

nm <- c(names(cdna), names(ncrna))
info <- data.frame(t(sapply(nm, function(w) {
  w <- strsplit(w, " ")[[1]]
  transcript <- w[1]
  gene <- gsub("gene:", "", w[grep("^gene:", w)])
  position <- gsub("chromosome:", "", w[grep("^chromosome:", w)])
  if (length(position) == 0) position <- gsub("scaffold:", "", w[grep("^scaffold:", w)])
  symbol <- gsub("gene_symbol:", "", w[grep("^gene_symbol", w)])
  c(transcript = transcript, 
    gene = gene,
    genome = strsplit(position, ":")[[1]][1], 
    chromosome = strsplit(position, ":")[[1]][2],
    start = strsplit(position, ":")[[1]][3],
    end = strsplit(position, ":")[[1]][4],
    strand = ifelse(strsplit(position, ":")[[1]][5] == "1", "+", "-"),
    symbol = symbol
  )
})), stringsAsFactors = FALSE)
rownames(info) <- NULL
info$start <- as.numeric(info$start)
info$end <- as.numeric(info$end)

txgr <- GRanges(seqnames = info$chromosome, 
                ranges = IRanges(start = info$start, end = info$end), 
                strand = info$strand)
mcols(txgr) <- info[, c("transcript", "gene", "genome", "symbol")]

erccgtf <- import(ercc_gtf, format = "gtf")
ercctx <- erccgtf
mcols(ercctx) <- data.frame(transcript = ercctx$gene_id,
                            gene = ercctx$gene_id,
                            genome = ercctx$source,
                            symbol = ercctx$gene_id)

txgr <- suppressWarnings(c(txgr, ercctx))
names(txgr) <- txgr$transcript

geneinfo <- info %>% group_by(gene) %>% 
  summarise(genome = unique(genome),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            strand = unique(strand),
            symbol = unique(symbol)) %>% 
  ungroup()

ggr <- GRanges(seqnames = geneinfo$chromosome,
               ranges = IRanges(start = geneinfo$start, end = geneinfo$end),
               strand = geneinfo$strand)
mcols(ggr) <- geneinfo[, c("gene", "genome", "symbol")]

erccgene <- erccgtf
mcols(erccgene) <- data.frame(gene = erccgene$gene_id,
                              genome = erccgene$source,
                              symbol = erccgene$gene_id)
ggr <- suppressWarnings(c(ggr, erccgene))
names(ggr) <- ggr$gene

gene_granges <- ggr
tx_granges <- txgr
saveRDS(list(gene_granges = gene_granges, tx_granges = tx_granges),
        file = gsub("cdna.all.fa$", paste0(mouse_ensembl_version, 
                                           ".cdna.ncrna.ercc92.granges.rds"), mouse_cdna_fa))

## Zebrafish
cdna <- readDNAStringSet(zebrafish_cdna_fa)
ncrna <- readDNAStringSet(zebrafish_ncrna_fa)

nm <- c(names(cdna), names(ncrna))
info <- data.frame(t(sapply(nm, function(w) {
  w <- strsplit(w, " ")[[1]]
  transcript <- w[1]
  gene <- gsub("gene:", "", w[grep("^gene:", w)])
  position <- gsub("chromosome:", "", w[grep("^chromosome:", w)])
  if (length(position) == 0) position <- gsub("scaffold:", "", w[grep("^scaffold:", w)])
  symbol <- gsub("gene_symbol:", "", w[grep("^gene_symbol", w)])
  c(transcript = transcript, 
    gene = gene,
    genome = strsplit(position, ":")[[1]][1], 
    chromosome = strsplit(position, ":")[[1]][2],
    start = strsplit(position, ":")[[1]][3],
    end = strsplit(position, ":")[[1]][4],
    strand = ifelse(strsplit(position, ":")[[1]][5] == "1", "+", "-"),
    symbol = symbol
  )
})), stringsAsFactors = FALSE)
rownames(info) <- NULL
info$start <- as.numeric(info$start)
info$end <- as.numeric(info$end)

txgr <- GRanges(seqnames = info$chromosome, 
                ranges = IRanges(start = info$start, end = info$end), 
                strand = info$strand)
mcols(txgr) <- info[, c("transcript", "gene", "genome", "symbol")]

erccgtf <- import(ercc_gtf, format = "gtf")
ercctx <- erccgtf
mcols(ercctx) <- data.frame(transcript = ercctx$gene_id,
                            gene = ercctx$gene_id,
                            genome = ercctx$source,
                            symbol = ercctx$gene_id)

txgr <- suppressWarnings(c(txgr, ercctx))
names(txgr) <- txgr$transcript

geneinfo <- info %>% group_by(gene) %>% 
  summarise(genome = unique(genome),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            strand = unique(strand),
            symbol = unique(symbol)) %>% 
  ungroup()

ggr <- GRanges(seqnames = geneinfo$chromosome,
               ranges = IRanges(start = geneinfo$start, end = geneinfo$end),
               strand = geneinfo$strand)
mcols(ggr) <- geneinfo[, c("gene", "genome", "symbol")]

erccgene <- erccgtf
mcols(erccgene) <- data.frame(gene = erccgene$gene_id,
                              genome = erccgene$source,
                              symbol = erccgene$gene_id)
ggr <- suppressWarnings(c(ggr, erccgene))
names(ggr) <- ggr$gene

gene_granges <- ggr
tx_granges <- txgr
saveRDS(list(gene_granges = gene_granges, tx_granges = tx_granges),
        file = gsub("cdna.all.fa$", paste0(zebrafish_ensembl_version, 
                                           ".cdna.ncrna.ercc92.granges.rds"), zebrafish_cdna_fa))

