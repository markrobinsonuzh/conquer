# conquer

This repository contains the code used to prepare the *conquer* (**con**sistent **qu**antification of **e**xternal **R**NA-seq data sets) repository, which provides access to consistently processed, analysis-ready single-cell RNA-seq data sets, together with quality control and exploratory analysis reports to help the user determine whether a particular data set fits their purposes. 

The repository is available at [http://imlspenticton.uzh.ch:3838/conquer/](http://imlspenticton.uzh.ch:3838/conquer/). A more detailed description of the workflow that was used to process the data is available in the publication listed below, within the application (see the "About" tab), and [here](shiny-download/about_conquer.md).

If you use *conquer*, please cite

* C Soneson & MD Robinson: [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612). Nature Methods (2018).

## Adapters

In cases where adapter trimming was necessary, we used a fasta filew ith adapters obtained from [https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt](https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt) (downloaded on July 25, 2016).

## Contact

If you have questions or comments about *conquer*, please contact [Charlotte Soneson](mailto:charlotte.soneson@uzh.ch) or [Mark D Robinson](mailto:mark.robinson@imls.uzh.ch).
