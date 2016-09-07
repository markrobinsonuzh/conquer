source("02_process_data.R")

salmonbin <- "/usr/local/software/SalmonBeta-0.6.1_DebianSqueeze/bin/salmon"
fastqcbin <- "FastQC_v0.11.6.devel/fastqc"
cutadaptbin <- "cutadapt"
Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap <- 
  readRDS("reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap.rds")
Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap <- 
  readRDS("reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap.rds")
Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges <- 
  readRDS("reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges.rds")
Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges <- 
  readRDS("reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges.rds")

adapterseq <- "AGATCGGAAGAGC"
#adapterseq <- "file:adapters.fa"

## E-MTAB-2805
process_data(id = "EMTAB2805", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "LibraryName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = FALSE, phenofile = "data-raw/EMTAB2805/EMTAB2805_pheno.txt",
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "cell_cycle_stage", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 25599176, datalink = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805/")

## GSE77847
process_data(id = "GSE77847", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "characteristics_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27016502, datalink = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77847")

## GSE44183human
process_data(id = "GSE44183-GPL11154", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183")

## GSE44183mouse
process_data(id = "GSE44183-GPL13112", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183")

## GSE44183mouse, trimmed
process_data(id = "GSE44183-GPL13112-trimmed", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = TRUE, adapterseq = "AGATCGGAAGAGC", cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183")

## GSE45719
process_data(id = "GSE45719", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "U", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 24408435, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719")

## GSE74596
process_data(id = "GSE74596", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
             libtype = "U", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27089380, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74596")

## GSE63818, HiSeq2500
process_data(id = "GSE63818-GPL16791", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("source_name_ch1", "characteristics_ch1"), 
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 26046443, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63818")

## GSE60749, HiSeq2000
process_data(id = "GSE60749-GPL13112", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.k15.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("source_name_ch1", "characteristics_ch1.1"), 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 25471879, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60749")

## GSE60749, HiSeq2500
process_data(id = "GSE60749-GPL17021", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.k15.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("characteristics_ch1", "characteristics_ch1.1"), 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 25471879, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60749")

## GSE57872
process_data(id = "GSE57872", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.k15.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("characteristics_ch1", "characteristics_ch1.1"), 
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24925914, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872")
