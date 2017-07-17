source("02_process_data.R")

#salmonbin <- "/usr/local/software/SalmonBeta-0.6.1_DebianSqueeze/bin/salmon"
#salmonbin <- "software/Salmon-0.7.2_linux_x86_64/bin/salmon"
salmonbin <- "/home/charlotte/software/Salmon-0.8.2_linux_x86_64/bin/salmon"
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

## E-MTAB-3929
process_data(id = "EMTAB3929", rtype = "single",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.8.2.sidx",
             libtype = "A", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "LibraryName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = FALSE, phenofile = "data-raw/EMTAB3929/EMTAB3929_pheno.txt",
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("Characteristics.developmental.stage.", "Characteristics.treatment."),
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27062923, datalink = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3929/",
             shortname = "Petropoulos2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".",
             bias = TRUE)

## GSE66507
process_data(id = "GSE66507", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.7.2.sidx",
             libtype = "A", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "characteristics_ch1.1", organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 26293300, datalink = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66507",
             shortname = "Blakeley2015", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE71585, HiSeq2000
process_data(id = "GSE71585-GPL13112", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.7.2.sidx",
             libtype = "A", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("characteristics_ch1.1"), 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 26727548, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",
             shortname = "Tasic2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE71585, HiSeq2500
process_data(id = "GSE71585-GPL17021", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.7.2.sidx",
             libtype = "A", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("characteristics_ch1.1"), 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 26727548, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",
             shortname = "Tasic2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE66688
process_data(id = "GSE66688", rtype = "paired", 
             index = "reference-files/danio-rerio/Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.0.7.2.sidx", 
             libtype = "A", salmonbin = "software/Salmon-0.7.2_linux_x86_64/bin/salmon",
             fastqcbin = fastqcbin, sncol = "SampleName", 
             txgenemap = Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Danio_rerio.GRCz10.87.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "characteristics_ch1", organism = "Danio rerio", genome = "GRCz10.87",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 25867923, datalink = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66688",
             shortname = "Satija2015", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## E-MTAB-2805
process_data(id = "EMTAB2805", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "LibraryName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = FALSE, phenofile = "data-raw/EMTAB2805/EMTAB2805_pheno.txt",
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "cell_cycle_stage", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 25599176, datalink = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805/",
             shortname = "Buettner2015", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE77847
process_data(id = "GSE77847", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "characteristics_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27016502, datalink = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77847",
             shortname = "Meyer2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE44183human
process_data(id = "GSE44183-GPL11154", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183",
             shortname = "Xue2013", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE44183mouse
process_data(id = "GSE44183-GPL13112", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183",
             shortname = "Xue2013", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE44183mouse, trimmed
process_data(id = "GSE44183-GPL13112-trimmed", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "Sample_Name_s",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = TRUE, adapterseq = "AGATCGGAAGAGC", cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 23892778, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183",
             shortname = "Xue2013", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE45719
process_data(id = "GSE45719", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "U", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 24408435, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719",
             shortname = "Deng2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE74596
process_data(id = "GSE74596", rtype = "single",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "U", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 27089380, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74596",
             shortname = "Engel2016", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE63818, HiSeq2500
process_data(id = "GSE63818-GPL16791", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = c("source_name_ch1", "characteristics_ch1"), 
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             pmid = 26046443, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63818",
             shortname = "Guo2015", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE60749, HiSeq2000 (25bp reads)
process_data(id = "GSE60749-GPL13112", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.k15.0.6.1.sidx",
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
             pmid = 25471879, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60749",
             shortname = "Kumar2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE60749, HiSeq2500 (25bp reads)
process_data(id = "GSE60749-GPL17021", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.k15.0.6.1.sidx",
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
             pmid = 25471879, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60749",
             shortname = "Kumar2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE57872 (25bp reads)
process_data(id = "GSE57872", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.k15.0.6.1.sidx",
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
             pmid = 24925914, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872",
             shortname = "Patel2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE41265
process_data(id = "GSE41265", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 1, lps = "bottom",
             pmid = 23685454, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41265",
             shortname = "Shalek2013", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE52529-GPL11154
process_data(id = "GSE52529-GPL11154", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24658644, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52529",
             shortname = "Trapnell2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE52529-GPL16791
process_data(id = "GSE52529-GPL16791", rtype = "paired",
             index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Homo sapiens", genome = "GRCh38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24658644, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52529",
             shortname = "Trapnell2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE48968-GPL13112
process_data(id = "GSE48968-GPL13112", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24919153, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48968",
             shortname = "Shalek2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE48968-GPL17021-125bp
process_data(id = "GSE48968-GPL17021-125bp", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24919153, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48968",
             shortname = "Shalek2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)

## GSE48968-GPL17021-25bp
process_data(id = "GSE48968-GPL17021-25bp", rtype = "paired",
             index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.k15.0.6.1.sidx",
             libtype = "IU", salmonbin = salmonbin,
             fastqcbin = fastqcbin, sncol = "SampleName",
             txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
             geodata = TRUE, phenofile = NULL,
             gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
             tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
             groupid = "source_name_ch1", 
             organism = "Mus musculus", genome = "GRCm38.84",
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp",
             nrw = 3, lps = "bottom",
             pmid = 24919153, datalink = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48968",
             shortname = "Shalek2014", multiqcbin = multiqcbin, 
             aspects = c("fastqc", "salmon", "multiqc", "mae", "scater"), topdir = ".", 
             bias = FALSE)
