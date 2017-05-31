source("02_process_data.R")

#salmonbin <- "/usr/local/software/SalmonBeta-0.6.1_DebianSqueeze/bin/salmon"
salmonbin <- "software/Salmon-0.7.2_linux_x86_64/bin/salmon"
fastqcbin <- "FastQC_v0.11.6.devel/fastqc"
cutadaptbin <- "cutadapt"
multiqcbin <- "/home/charlotte/miniconda2/bin/multiqc"

Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap <- 
  readRDS("reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap.rds")
Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap <- 
  readRDS("reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap.rds")
Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.txgenemap <- 
  readRDS("reference-files/danio-rerio/Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.txgenemap.rds")

Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges <- 
  readRDS("reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds")
Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges <- 
  readRDS("reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges.rds")
Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.granges <- 
  readRDS("reference-files/danio-rerio/Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.granges.rds")

adapterseq <- "AGATCGGAAGAGC"
#adapterseq <- "file:adapters.fa"

## GSE84527
process_data(id = "GSE84527", rtype = "paired", 
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.7.2.sidx", 
             libtype = "A", salmonbin = "software/Salmon-0.7.2_linux_x86_64/bin/salmon",
             fastqcbin = fastqcbin, sncol = "SampleName", 
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "characteristics_ch1.2", organism = "Mus Musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27716482, datalink = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84527",
             shortname = "Aihara2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"),
             topdir = "bulk_rnaseq")

## PRJNA281360
process_data(id = "PRJNA281360", rtype = "paired", 
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.7.2.sidx", 
             libtype = "A", salmonbin = "software/Salmon-0.7.2_linux_x86_64/bin/salmon",
             fastqcbin = fastqcbin, sncol = "Experiment", 
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = FALSE, phenofile = "bulk_rnaseq/data-raw/PRJNA281360/PRJNA281360_pheno.txt",
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("CellType", "Sex"), organism = "Mus Musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = NA, datalink = "https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=281360",
             shortname = "ImmGen", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"),
             topdir = "bulk_rnaseq")

## E-GEUV-1
process_data(id = "EGEUV1", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.7.2.sidx",
             libtype = "A", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Experiment",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = FALSE, phenofile = "bulk_rnaseq/data-raw/EGEUV1/EGEUV1_pheno.txt",
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("center", "population"), organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = NA, datalink = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/",
             shortname = "Geuvadis", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"),
             topdir = "bulk_rnaseq")
