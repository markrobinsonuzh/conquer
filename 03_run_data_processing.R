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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

# # GSE29087 (protocol designed to capture 5' ends of transcripts)
# process_data(id = "GSE29087", rtype = "single",
#              index = "reference-files/homo-sapiens/Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.sidx",
#              libtype = "U", salmonbin = salmonbin,
#              fastqcbin = fastqcbin, sncol = "SampleName",
#              txgenemap = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.txgenemap,
#              geodata = TRUE, phenofile = NULL,
#              gene_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$gene_granges,
#              tx_granges = Homo_sapiens.GRCh38.84.cdna.ncrna.ercc92.granges$tx_granges,
#              groupid = "source_name_ch1", organism = "Homo sapiens", genome = "GRCh38.84",
#              dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

## GSE41265 (no comparison to be made, only one group)
# process_data(id = "GSE41265", rtype = "paired",
#              index = "reference-files/mus-musculus/Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.sidx",
#              libtype = "IU", salmonbin = salmonbin,
#              fastqcbin = fastqcbin, sncol = "SampleName",
#              txgenemap = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.txgenemap,
#              geodata = TRUE, phenofile = NULL,
#              gene_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$gene_granges,
#              tx_granges = Mus_musculus.GRCm38.84.cdna.ncrna.ercc92.granges$tx_granges,
#              groupid = "characteristics_ch1.1", organism = "Mus musculus", genome = "GRCm38.84",
#              dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = TRUE, adapterseq = "AGATCGGAAGAGC", cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")

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
             dotrim = FALSE, adapterseq = NULL, cutadaptbin = cutadaptbin, tmp_dir = "tmp")
