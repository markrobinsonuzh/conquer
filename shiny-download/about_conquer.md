# About *conquer*

The *conquer* (**con**sistent **qu**antification of **e**xternal **r**na-seq data) repository is developed by [Charlotte Soneson](mailto:charlotte.soneson@imls.uzh.ch) and [Mark D Robinson](mailto:mark.robinson@imls.uzh.ch) at the University of Zurich, Switzerland. It is implemented in [shiny](http://shiny.rstudio.com/) and provides access to consistently processed public single-cell RNA-seq data sets. Below is a short description of the workflow used to process the raw reads in order to generate the data provided in the repository.

<a name="indexing"/>
</a>
## Index building

In order to use [Salmon](https://combine-lab.github.io/salmon/) to quantify the transcript abundances in a given sample, we first need to index the corresponding reference transcriptome. For a given organism, we download the fasta files containing cDNA and ncRNA sequences from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html), complement these with [ERCC spike-in sequences](https://www.thermofisher.com/order/catalog/product/4456739), and build a Salmon  quasi-mapping index for the entire catalog. Note that the [scater](https://bioconductor.org/packages/devel/bioc/html/scater.html) report for a given data set (available in the **scater_html** column) details the precise version of the transcriptome that was used for the quantification. For data sets with "long" reads (longer than 50 bp) we use the default `k=31`, while for "short reads" (typically around 25 bp) we set `k=15`. 

We also create a lookup table relating transcript IDs to the corresponding gene IDs. This information is obtained by parsing the sequence names in the cDNA and ncRNA fasta files. From these names we also obtain the genomic coordinates for each feature. 

## Sample list and run matching

The first step is to determine the set of samples included in a given data set. We download a "RunInfo.csv" file for the data set from [SRA](http://www.ncbi.nlm.nih.gov/sra) and a Series Matrix file from [GEO](http://www.ncbi.nlm.nih.gov/geo/), in order to link samples both to individual runs and to phenotypic information. If the data set is not available from GEO, we construct a phenotype data file from the information provided by the corresponding repository. 

## Quality control

For each sample in the data set, we find all the corresponding runs, and download and concatenate the corresponding FastQ files from SRA. There is also an optional step to trim adapters from the reads using [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html). Next, we run [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate a quality control file for each concatenated read file (one or two files per sample depending on whether it was processed with a single-end or paired-end sequencing protocol). 

## Abundance quantification

After the QC, we run [Salmon](https://combine-lab.github.io/salmon/) to estimate the abundance of each transcript from the catalog described [above](#indexing) in each sample. The Salmon output files are then compressed in an archive and can be downloaded from *conquer* (see the **salmon_tar** column).

## Summary report - MultiQC

Once FastQC and Salmon have been applied to all samples in the data set, we run [MultiQC](http://multiqc.info/) to summarise all the information into one report. This can also be downloaded from *conquer* (see the **MultiQC_html** column). This report contains quality scores for all the samples and can be used to determine if there are problematic samples and whether the data set is good enough for the purposes of the user or needs to be subsetted.

## Data summarisation

The abundances estimated by Salmon are summarised and provided to the user via *conquer* in the form of a [MultiAssayExperiment](https://bioconductor.org/packages/devel/bioc/html/MultiAssayExperiment.html) object. This object can be downloaded via the buttons in the **MultiAssayExperiment** column. To generate this object, we first use the [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) package to read the Salmon output into R. This returns both count estimates and TPM estimates for each transcript. Next, we summarise the transcript-level information to the gene level. The gene-level TPM is defined as the sum of the TPMs of the corresponding transcripts, and similarly for the gene-level counts. We also provide "scaled TPMs" (see [http://f1000research.com/articles/4-1521/](http://f1000research.com/articles/4-1521/) or the [tximport vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) for a discussion), that is, summarised TPMs scaled to a "count scale". In the summarisation step, we make use of the transcript-to-gene lookup table generated [above](#indexing).

The provided MultiAssayExperiment object contains two "experiments", corresponding to the gene-level and transcript-level values. The gene-level experiment contains four "assays": 

- *TPM*
- *count*
- *count_lstpm* (count-scale length-scaled TPMs)
- *avetxlength* (the average transcript length, which can be used as offsets in count models based on the *count* assay, see [http://f1000research.com/articles/4-1521/](http://f1000research.com/articles/4-1521/)).  

The transcript-level experiment contains three "assays":

- *TPM*
- *count*
- *efflength* (the effective length estimated by Salmon)

The MultiAssayExperiment also contains the phenotypic data (in the *pData* slot), as well as some metadata for the data set (the genome, the organism and the Salmon index that was used for the quantification). 

## Summary report - scater

In order to give users another way of investigating whether a data set is useful for their purposes, we also provide an exploratory analysis report. This is largely based on functions from the [scater](https://bioconductor.org/packages/devel/bioc/html/scater.html) Bioconductor package, applied to data extracted from the MultiAssayExperiment object. The report calculates and visualises various quality measures for the cells, and provides low-dimensional representations of the cells, colored by different phenotypic annotations. 

## Acknowledgements

We would like to thank Simon Andrews for help with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), Mike Love and Valentine Svensson for providing instructions for how to retrieve the URL for the FastQ file(s) of a given SRA run (see [here](https://gist.github.com/mikelove/f539631f9e187a8931d34779436a1c01) and [here](http://www.nxn.se/valent/streaming-rna-seq-data-from-ena)), Davis McCarthy for input regarding [scater](https://bioconductor.org/packages/devel/bioc/html/scater.html) and Nicholas Hamilton for instructions on how to generate a standardized report based on a provided R object (see [here](http://stackoverflow.com/questions/37097535/generate-report-in-r))). Finally, we would like to acknowledge the developers of all the tools we use to prepare the data for *conquer*. 

## Presentations/publications

*conquer* was presented as a [poster](http://imlspenticton.uzh.ch/robinson_lab/conquer/presentation/Soneson-poster-singlecellgenomics-2016.pdf) at the [Single Cell Genomics](https://coursesandconferences.wellcomegenomecampus.org/events/item.aspx?e=596) conference in Hinxton, UK, in September 2016.

## Code

The code used for *conquer* is available via [GitHub](https://github.com/markrobinsonuzh/conquer).