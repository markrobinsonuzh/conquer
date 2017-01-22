# Using the conquer database
  
  

To use a data set provided in the *conquer* database, download the corresponding
**R** object from the **MultiAssayExperiment** column. As an illustration, we 
will assume that the file for the *GSE41265* data set has been downloaded and is 
available in the current working directory. First, load the 
[SummarizedExperiment](http://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
and
[MultiAssayExperiment](http://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)
packages and read the file into R:


```r
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
(gse41265 <- readRDS("GSE41265.rds"))
```

```
## A MultiAssayExperiment object of 2 listed
##  experiments with user-defined names and respective classes. 
##  Containing an ExperimentList class object of length 2: 
##  [1] gene: RangedSummarizedExperiment with 45686 rows and 18 columns 
##  [2] tx: RangedSummarizedExperiment with 113560 rows and 18 columns 
## To access: 
##  experiments() - to obtain the ExperimentList instance 
##  pData() - for the primary/phenotype DataFrame 
##  sampleMap() - for the sample availability DataFrame 
##  metadata() - for the metadata object of ANY class 
## See also: subsetByAssay(), subsetByRow(), subsetByColumn()
```

The resulting object contains both gene and transcript abundances. 


```r
experiments(gse41265)
```

```
## ExperimentList class object of length 2: 
##  [1] gene: RangedSummarizedExperiment with 45686 rows and 18 columns 
##  [2] tx: RangedSummarizedExperiment with 113560 rows and 18 columns
```

## Gene-level data

To access the gene abundances, get the **gene** experiment:


```r
(gse41265_gene <- experiments(gse41265)[["gene"]])
```

```
## class: RangedSummarizedExperiment 
## dim: 45686 18 
## metadata(0):
## assays(4): TPM count count_lstpm avetxlength
## rownames(45686): ENSMUSG00000000001.4 ENSMUSG00000000003.15 ...
##   ERCC-00170 ERCC-00171
## rowData names(3): gene genome symbol
## colnames(18): GSM1012777 GSM1012778 ... GSM1012793 GSM1012794
## colData names(0):
```

This object contains four slots, which can be accessed via the **assays** function:

- **TPM**: transcripts per million abundance estimates for each gene, obtained
by summing the transcript TPMs for the gene's isoforms.
- **count**: gene read counts, obtained by summing the estimated read counts for
the gene's isoforms.
- **count_lstpm**: length-scaled TPMs, which provide and alternative abundance
measure on the "count scale", which is not correlated with the average
transcript length in a given sample. See the
[tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html)
Bioconductor package for more information.
- **avetxlength**: the average length of the transcripts expressed in each
sample for each gene. See the
[tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html)
Bioconductor package for more information.

Each of these slots is a matrix with the respective values for each gene and
each sample.


```r
head(assays(gse41265_gene)[["TPM"]])
```

```
##                       GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUSG00000000001.4     16.9260 68.9189000 27.2136000    25.6035
## ENSMUSG00000000003.15     0.0000  0.0000000  0.0000000     0.0000
## ENSMUSG00000000028.14    84.6329  0.0354825  0.0261671     0.0000
## ENSMUSG00000000031.15     0.0000  0.0000000  0.0000000     0.0000
## ENSMUSG00000000037.16    21.0225  0.0000000  0.0000000     0.0000
## ENSMUSG00000000049.11     0.0000  0.0000000  0.0000000     0.0000
##                       GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUSG00000000001.4    1.003940 13.9319000    1.62186   27.26950
## ENSMUSG00000000003.15   0.000000  0.0000000    0.00000    0.00000
## ENSMUSG00000000028.14   0.028031  0.0816208    0.00000   18.28858
## ENSMUSG00000000031.15   0.000000  0.0000000    0.00000    0.00000
## ENSMUSG00000000037.16   0.000000  0.0000000    0.00000    0.00000
## ENSMUSG00000000049.11   0.000000  0.0000000    0.00000    0.00000
##                       GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUSG00000000001.4    0.220791    0.14094    21.8415  42.991200
## ENSMUSG00000000003.15   0.000000    0.00000     0.0000   0.000000
## ENSMUSG00000000028.14  58.226400    0.00000     0.0000   0.118938
## ENSMUSG00000000031.15   0.000000    0.00000     0.0000   0.000000
## ENSMUSG00000000037.16   0.000000    0.00000     0.0000   0.000000
## ENSMUSG00000000049.11   0.000000    0.00000     0.0000   0.000000
##                       GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUSG00000000001.4     156.447   0.109962  1.3144600 27.3401000
## ENSMUSG00000000003.15      0.000   0.000000  0.0000000  0.0000000
## ENSMUSG00000000028.14      0.000   3.537801  0.0382929  0.0340262
## ENSMUSG00000000031.15      0.000   0.000000  0.0000000  0.0000000
## ENSMUSG00000000037.16      0.000   0.000000  0.0000000  0.5705780
## ENSMUSG00000000049.11      0.000   0.000000  0.0000000  0.8658010
##                       GSM1012793 GSM1012794
## ENSMUSG00000000001.4     23.0052 32.8455000
## ENSMUSG00000000003.15     0.0000  0.0000000
## ENSMUSG00000000028.14     0.0000  0.0453766
## ENSMUSG00000000031.15     0.0000  0.0000000
## ENSMUSG00000000037.16     0.0000  0.0000000
## ENSMUSG00000000049.11     0.0000  0.0000000
```

```r
head(assays(gse41265_gene)[["count"]])
```

```
##                       GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUSG00000000001.4     777.612    3903.94    2088.36    1493.55
## ENSMUSG00000000003.15      0.000       0.00       0.00       0.00
## ENSMUSG00000000028.14   2352.123       1.00       1.00       0.00
## ENSMUSG00000000031.15      0.000       0.00       0.00       0.00
## ENSMUSG00000000037.16     90.000       0.00       0.00       0.00
## ENSMUSG00000000049.11      0.000       0.00       0.00       0.00
##                       GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUSG00000000001.4     71.8391 875.125000    93.6992   1645.560
## ENSMUSG00000000003.15     0.0000   0.000000     0.0000      0.000
## ENSMUSG00000000028.14     1.0000   2.000001     0.0000    672.918
## ENSMUSG00000000031.15     0.0000   0.000000     0.0000      0.000
## ENSMUSG00000000037.16     0.0000   0.000000     0.0000      0.000
## ENSMUSG00000000049.11     0.0000   0.000000     0.0000      0.000
##                       GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUSG00000000001.4       11.00          9    1515.99    1843.77
## ENSMUSG00000000003.15       0.00          0       0.00       0.00
## ENSMUSG00000000028.14    1823.75          0       0.00       1.00
## ENSMUSG00000000031.15       0.00          0       0.00       0.00
## ENSMUSG00000000037.16       0.00          0       0.00       0.00
## ENSMUSG00000000049.11       0.00          0       0.00       0.00
##                       GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUSG00000000001.4     9756.93     6.0000    69.0064    1614.74
## ENSMUSG00000000003.15       0.00     0.0000     0.0000       0.00
## ENSMUSG00000000028.14       0.00   121.1216     1.0000       1.00
## ENSMUSG00000000031.15       0.00     0.0000     0.0000       0.00
## ENSMUSG00000000037.16       0.00     0.0000     0.0000      30.00
## ENSMUSG00000000049.11       0.00     0.0000     0.0000      16.00
##                       GSM1012793 GSM1012794
## ENSMUSG00000000001.4     1277.01     1459.1
## ENSMUSG00000000003.15       0.00        0.0
## ENSMUSG00000000028.14       0.00        1.0
## ENSMUSG00000000031.15       0.00        0.0
## ENSMUSG00000000037.16       0.00        0.0
## ENSMUSG00000000049.11       0.00        0.0
```

```r
head(assays(gse41265_gene)[["count_lstpm"]])
```

```
##                       GSM1012777   GSM1012778   GSM1012779 GSM1012780
## ENSMUSG00000000001.4    752.8798 3720.7822146 2032.4459094   1452.948
## ENSMUSG00000000003.15     0.0000    0.0000000    0.0000000      0.000
## ENSMUSG00000000028.14  1872.5288    0.9528578    0.9720902      0.000
## ENSMUSG00000000031.15     0.0000    0.0000000    0.0000000      0.000
## ENSMUSG00000000037.16   290.8873    0.0000000    0.0000000      0.000
## ENSMUSG00000000049.11     0.0000    0.0000000    0.0000000      0.000
##                       GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUSG00000000001.4  69.9132476 847.647213   91.01965  1582.7900
## ENSMUSG00000000003.15  0.0000000   0.000000    0.00000     0.0000
## ENSMUSG00000000028.14  0.9709756   2.470152    0.00000   528.0124
## ENSMUSG00000000031.15  0.0000000   0.000000    0.00000     0.0000
## ENSMUSG00000000037.16  0.0000000   0.000000    0.00000     0.0000
## ENSMUSG00000000049.11  0.0000000   0.000000    0.00000     0.0000
##                       GSM1012785 GSM1012786 GSM1012787  GSM1012788
## ENSMUSG00000000001.4    10.61746   8.385097   1487.373 1497.153135
## ENSMUSG00000000003.15    0.00000   0.000000      0.000    0.000000
## ENSMUSG00000000028.14 1392.76246   0.000000      0.000    2.060275
## ENSMUSG00000000031.15    0.00000   0.000000      0.000    0.000000
## ENSMUSG00000000037.16    0.00000   0.000000      0.000    0.000000
## ENSMUSG00000000049.11    0.00000   0.000000      0.000    0.000000
##                       GSM1012789 GSM1012790 GSM1012791   GSM1012792
## ENSMUSG00000000001.4    9318.475   5.840962 66.3240941 1559.5448194
## ENSMUSG00000000003.15      0.000   0.000000  0.0000000    0.0000000
## ENSMUSG00000000028.14      0.000  93.474485  0.9610815    0.9654487
## ENSMUSG00000000031.15      0.000   0.000000  0.0000000    0.0000000
## ENSMUSG00000000037.16      0.000   0.000000  0.0000000   10.1246935
## ENSMUSG00000000049.11      0.000   0.000000  0.0000000   15.4503225
##                       GSM1012793   GSM1012794
## ENSMUSG00000000001.4    1201.013 1407.2935134
## ENSMUSG00000000003.15      0.000    0.0000000
## ENSMUSG00000000028.14      0.000    0.9670719
## ENSMUSG00000000031.15      0.000    0.0000000
## ENSMUSG00000000037.16      0.000    0.0000000
## ENSMUSG00000000049.11      0.000    0.0000000
```

```r
head(assays(gse41265_gene)[["avetxlength"]])
```

```
##                       GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUSG00000000001.4   3026.4100  3015.1200  3017.9000  3018.6900
## ENSMUSG00000000003.15   553.6516   553.6516   553.6516   553.6516
## ENSMUSG00000000028.14  1830.7909  1500.1200  1502.9000  1462.8034
## ENSMUSG00000000031.15  1022.8460  1022.8460  1022.8460  1022.8460
## ENSMUSG00000000037.16   282.0180   870.1116   870.1116   870.1116
## ENSMUSG00000000049.11   943.5560   943.5560   943.5560   943.5560
##                       GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUSG00000000001.4   3021.2400  3018.8100  3015.9000  3014.2200
## ENSMUSG00000000003.15   553.6516   553.6516   553.6516   553.6516
## ENSMUSG00000000028.14  1506.2400  1177.6154  1462.8034  1837.8923
## ENSMUSG00000000031.15  1022.8460  1022.8460  1022.8460  1022.8460
## ENSMUSG00000000037.16   870.1116   870.1116   870.1116   870.1116
## ENSMUSG00000000049.11   943.5560   943.5560   943.5560   943.5560
##                       GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUSG00000000001.4   3013.6100  3005.7400  3017.4300  3022.5500
## ENSMUSG00000000003.15   553.6516   553.6516   553.6516   553.6516
## ENSMUSG00000000028.14  1894.6100  1462.8034  1462.8034   592.5530
## ENSMUSG00000000031.15  1022.8460  1022.8460  1022.8460  1022.8460
## ENSMUSG00000000037.16   870.1116   870.1116   870.1116   870.1116
## ENSMUSG00000000049.11   943.5560   943.5560   943.5560   943.5560
##                       GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUSG00000000001.4   3015.2300  3017.5800  3014.5600  3015.5600
## ENSMUSG00000000003.15   553.6516   553.6516   553.6516   553.6516
## ENSMUSG00000000028.14  1462.8034  1893.3784  1499.5600  1500.5600
## ENSMUSG00000000031.15  1022.8460  1022.8460  1022.8460  1022.8460
## ENSMUSG00000000037.16   870.1116   870.1116   870.1116  2684.5600
## ENSMUSG00000000049.11   943.5560   943.5560   943.5560   943.5560
##                       GSM1012793 GSM1012794
## ENSMUSG00000000001.4   3012.8000  3006.4800
## ENSMUSG00000000003.15   553.6516   553.6516
## ENSMUSG00000000028.14  1462.8034  1491.4800
## ENSMUSG00000000031.15  1022.8460  1022.8460
## ENSMUSG00000000037.16   870.1116   870.1116
## ENSMUSG00000000049.11   943.5560   943.5560
```

## Transcript-level data

To access the transcript abundances, get instead the **transcript** experiment:


```r
(gse41265_tx <- experiments(gse41265)[["tx"]])
```

```
## class: RangedSummarizedExperiment 
## dim: 113560 18 
## metadata(0):
## assays(3): TPM count efflength
## rownames(113560): ENSMUST00000178537.1 ENSMUST00000178862.1 ...
##   ERCC-00170 ERCC-00171
## rowData names(4): transcript gene genome symbol
## colnames(18): GSM1012777 GSM1012778 ... GSM1012793 GSM1012794
## colData names(0):
```

This object contains three slots, which can be accessed via the **assays**
function:

- **TPM**: transcripts per million abundance estimates for each transcript.
- **count**: transcript read counts.
- **efflength**: effective transcript lengths.

Each of these slots is a matrix with the respective values for each transcript
and each sample.


```r
head(assays(gse41265_tx)[["TPM"]])
```

```
##                      GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012793 GSM1012794
## ENSMUST00000178537.1          0          0
## ENSMUST00000178862.1          0          0
## ENSMUST00000177564.1          0          0
## ENSMUST00000196221.1          0          0
## ENSMUST00000179664.1          0          0
## ENSMUST00000179520.1          0          0
```

```r
head(assays(gse41265_tx)[["count"]])
```

```
##                      GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUST00000178537.1          0          0          0          0
## ENSMUST00000178862.1          0          0          0          0
## ENSMUST00000177564.1          0          0          0          0
## ENSMUST00000196221.1          0          0          0          0
## ENSMUST00000179664.1          0          0          0          0
## ENSMUST00000179520.1          0          0          0          0
##                      GSM1012793 GSM1012794
## ENSMUST00000178537.1          0          0
## ENSMUST00000178862.1          0          0
## ENSMUST00000177564.1          0          0
## ENSMUST00000196221.1          0          0
## ENSMUST00000179664.1          0          0
## ENSMUST00000179520.1          0          0
```

```r
head(assays(gse41265_tx)[["efflength"]])
```

```
##                      GSM1012777 GSM1012778 GSM1012779 GSM1012780
## ENSMUST00000178537.1         12         12         12         12
## ENSMUST00000178862.1         14         14         14         14
## ENSMUST00000177564.1         16         16         16         16
## ENSMUST00000196221.1          9          9          9          9
## ENSMUST00000179664.1         11         11         11         11
## ENSMUST00000179520.1         11         11         11         11
##                      GSM1012781 GSM1012782 GSM1012783 GSM1012784
## ENSMUST00000178537.1         12         12         12         12
## ENSMUST00000178862.1         14         14         14         14
## ENSMUST00000177564.1         16         16         16         16
## ENSMUST00000196221.1          9          9          9          9
## ENSMUST00000179664.1         11         11         11         11
## ENSMUST00000179520.1         11         11         11         11
##                      GSM1012785 GSM1012786 GSM1012787 GSM1012788
## ENSMUST00000178537.1         12         12         12         12
## ENSMUST00000178862.1         14         14         14         14
## ENSMUST00000177564.1         16         16         16         16
## ENSMUST00000196221.1          9          9          9          9
## ENSMUST00000179664.1         11         11         11         11
## ENSMUST00000179520.1         11         11         11         11
##                      GSM1012789 GSM1012790 GSM1012791 GSM1012792
## ENSMUST00000178537.1         12         12         12         12
## ENSMUST00000178862.1         14         14         14         14
## ENSMUST00000177564.1         16         16         16         16
## ENSMUST00000196221.1          9          9          9          9
## ENSMUST00000179664.1         11         11         11         11
## ENSMUST00000179520.1         11         11         11         11
##                      GSM1012793 GSM1012794
## ENSMUST00000178537.1         12         12
## ENSMUST00000178862.1         14         14
## ENSMUST00000177564.1         16         16
## ENSMUST00000196221.1          9          9
## ENSMUST00000179664.1         11         11
## ENSMUST00000179520.1         11         11
```

## Sample annotations

The sample annotations, downloaded from GEO, are also available in the object:


```r
pdata <- Biobase::pData(gse41265)
head(pdata, 2)
```

```
## DataFrame with 2 rows and 47 columns
##                     title geo_accession                status
##                  <factor>      <factor>              <factor>
## GSM1012777 Single cell S1    GSM1012777 Public on May 19 2013
## GSM1012778 Single cell S2    GSM1012778 Public on May 19 2013
##            submission_date last_update_date     type channel_count
##                   <factor>         <factor> <factor>      <factor>
## GSM1012777     Oct 01 2012      May 19 2013      SRA             1
## GSM1012778     Oct 01 2012      May 19 2013      SRA             1
##               source_name_ch1 organism_ch1 characteristics_ch1
##                      <factor>     <factor>            <factor>
## GSM1012777 BMDC (4h LPS stim) Mus musculus     strain: C57BL/6
## GSM1012778 BMDC (4h LPS stim) Mus musculus     strain: C57BL/6
##                                           characteristics_ch1.1
##                                                        <factor>
## GSM1012777 cell type: Bone Marrow-derived Dendritic Cell (BMDC)
## GSM1012778 cell type: Bone Marrow-derived Dendritic Cell (BMDC)
##                 characteristics_ch1.2 characteristics_ch1.3
##                              <factor>              <factor>
## GSM1012777 treatment: LPS-stimulation    cell count: 1 cell
## GSM1012778 treatment: LPS-stimulation    cell count: 1 cell
##            characteristics_ch1.4
##                         <factor>
## GSM1012777                      
## GSM1012778                      
##                                                                               growth_protocol_ch1
##                                                                                          <factor>
## GSM1012777 Cells were cultured and stimulated with LPS as previously described (Amit et. al 2009)
## GSM1012778 Cells were cultured and stimulated with LPS as previously described (Amit et. al 2009)
##            molecule_ch1              extract_protocol_ch1
##                <factor>                          <factor>
## GSM1012777    polyA RNA cDNA synthesis and amplification:
## GSM1012778    polyA RNA cDNA synthesis and amplification:
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           extract_protocol_ch1.1
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         <factor>
## GSM1012777 We used the SMARTer Ultra Low RNA Kit (Clontech, Mountain View, CA) to prepare amplified cDNA. We added 1   l of 12   M 3' SMART primer (5'   AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N (N = A, C, G, or T; N-1 = A, G, or C)), 1   l of H2O, and 2.5   l of Reaction Buffer onto the RNA capture beads. We mixed them well by pipetting, heated the mixture at 72  C for 3 minutes and placed it on ice. First-strand cDNA was synthesized with this RNA primer mix by adding 2   l of 5x first-strand buffer, 0.25   l of 100mM DTT, 1   l of 10 mM dNTPs, 1   l of 12   M SMARTer II A Oligo (5'   AAGCAGTGGTATCAACGCAGAGTACXXXXX (X = undisclosed base in the proprietary SMARTer oligo sequence)), 100 U SMARTScribe RT, and 10 U RNase Inhibitor in a total volume of 10   l and incubating at 42  C for 90 minutes followed by 10 minutes at 70  C. We purified the first strand cDNA by adding 25   l of room temperature AMPure XP SPRI beads (Beckman Coulter Genomics, Danvers, MA), mixing well by pipetting, incubating at room temperature for 8 minutes. We removed the supernatant from the beads after a good separation was established. We carried out all of the above steps in a PCR product   free clean room. We amplified the cDNA by adding 5   l of 10x Advantage 2 PCR Buffer, 2   l of 10 mM dNTPs, 2   l of 12   M IS PCR primer (5'    AAGCAGTGGTATCAACGCAGAGT), 2   l of 50x Advantage 2 Polymerase Mix, and 39   l H2O in a total volume of 50   l. We performed the PCR at 95  C for 1 minute, followed by 21 cycles of 15 seconds at 95  C, 30 seconds at 65  C and 6 minutes at 68  C, followed by another 10 minutes at 72  C for final extension. We purified the amplified cDNA by adding 90   l of AMPure XP SPRI beads and washing with 80% ethanol.
## GSM1012778 We used the SMARTer Ultra Low RNA Kit (Clontech, Mountain View, CA) to prepare amplified cDNA. We added 1   l of 12   M 3' SMART primer (5'   AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N (N = A, C, G, or T; N-1 = A, G, or C)), 1   l of H2O, and 2.5   l of Reaction Buffer onto the RNA capture beads. We mixed them well by pipetting, heated the mixture at 72  C for 3 minutes and placed it on ice. First-strand cDNA was synthesized with this RNA primer mix by adding 2   l of 5x first-strand buffer, 0.25   l of 100mM DTT, 1   l of 10 mM dNTPs, 1   l of 12   M SMARTer II A Oligo (5'   AAGCAGTGGTATCAACGCAGAGTACXXXXX (X = undisclosed base in the proprietary SMARTer oligo sequence)), 100 U SMARTScribe RT, and 10 U RNase Inhibitor in a total volume of 10   l and incubating at 42  C for 90 minutes followed by 10 minutes at 70  C. We purified the first strand cDNA by adding 25   l of room temperature AMPure XP SPRI beads (Beckman Coulter Genomics, Danvers, MA), mixing well by pipetting, incubating at room temperature for 8 minutes. We removed the supernatant from the beads after a good separation was established. We carried out all of the above steps in a PCR product   free clean room. We amplified the cDNA by adding 5   l of 10x Advantage 2 PCR Buffer, 2   l of 10 mM dNTPs, 2   l of 12   M IS PCR primer (5'    AAGCAGTGGTATCAACGCAGAGT), 2   l of 50x Advantage 2 Polymerase Mix, and 39   l H2O in a total volume of 50   l. We performed the PCR at 95  C for 1 minute, followed by 21 cycles of 15 seconds at 95  C, 30 seconds at 65  C and 6 minutes at 68  C, followed by another 10 minutes at 72  C for final extension. We purified the amplified cDNA by adding 90   l of AMPure XP SPRI beads and washing with 80% ethanol.
##                                                                                 extract_protocol_ch1.2
##                                                                                               <factor>
## GSM1012777 We created Illumina sequencing libraries from this amplified cDNA using standard protocols.
## GSM1012778 We created Illumina sequencing libraries from this amplified cDNA using standard protocols.
##                             extract_protocol_ch1.3
##                                           <factor>
## GSM1012777 cDNA shearing and library construction:
## GSM1012778 cDNA shearing and library construction:
##                                                                                                                                                                                                                                                                                                                                         extract_protocol_ch1.4
##                                                                                                                                                                                                                                                                                                                                                       <factor>
## GSM1012777 We added the purification buffer (Clontech) to the amplified cDNA to make a total volume of 76   l. We sheared the cDNA in a 100   l tube with 10% Duty Cycle, 5% Intensity and 200 Cycles/Burst for 5 minutes in the frequency sweeping mode (Covaris S2 machine, Woburn, MA). We purified the sheared cDNA with 2.2 volumes AMPure XP SPRI beads.
## GSM1012778 We added the purification buffer (Clontech) to the amplified cDNA to make a total volume of 76   l. We sheared the cDNA in a 100   l tube with 10% Duty Cycle, 5% Intensity and 200 Cycles/Burst for 5 minutes in the frequency sweeping mode (Covaris S2 machine, Woburn, MA). We purified the sheared cDNA with 2.2 volumes AMPure XP SPRI beads.
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         extract_protocol_ch1.5
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       <factor>
## GSM1012777 We prepared indexed paired-end libraries for Illumina sequencing as described (J. Z. Levin et al., Nature Methods 7, 709 (2010)., with the following modifications. First, we used a different indexing adaptor (containing an 8-base barcode) for each library. Second, we size-selected the ligation product by using two rounds of 0.7 volume of AMPure XP SPRI bead cleanup with the first round starting volume at 100   l. Third, we performed PCR with Phusion High-Fidelity DNA polymerase with GC buffer and 2 M betaine. Fourth, we used 55  C as the annealing temperature in PCR with the universal indexing primers (forward primer 5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC, reverse primer 5'-CAAGCAGAAGACGGCATACGAGAT). Fifth, we performed 12 cycles of PCR. Sixth, we removed PCR primers using two rounds of 1.0 volume of AMPure beads.
## GSM1012778 We prepared indexed paired-end libraries for Illumina sequencing as described (J. Z. Levin et al., Nature Methods 7, 709 (2010)., with the following modifications. First, we used a different indexing adaptor (containing an 8-base barcode) for each library. Second, we size-selected the ligation product by using two rounds of 0.7 volume of AMPure XP SPRI bead cleanup with the first round starting volume at 100   l. Third, we performed PCR with Phusion High-Fidelity DNA polymerase with GC buffer and 2 M betaine. Fourth, we used 55  C as the annealing temperature in PCR with the universal indexing primers (forward primer 5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC, reverse primer 5'-CAAGCAGAAGACGGCATACGAGAT). Fifth, we performed 12 cycles of PCR. Sixth, we removed PCR primers using two rounds of 1.0 volume of AMPure beads.
##            taxid_ch1 description
##             <factor>    <factor>
## GSM1012777     10090          S1
## GSM1012778     10090          S2
##                                                                                                                                                                                                                                           data_processing
##                                                                                                                                                                                                                                                  <factor>
## GSM1012777 We created a Bowtie index based on the UCSC knownGene (8) transcriptome, and aligned paired-end reads directly to this index using Bowtie v 0.12.7 with command line options -q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 1000 -a -m 200.
## GSM1012778 We created a Bowtie index based on the UCSC knownGene (8) transcriptome, and aligned paired-end reads directly to this index using Bowtie v 0.12.7 with command line options -q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 1000 -a -m 200.
##                                                                                                                                                                                                                                       data_processing.1
##                                                                                                                                                                                                                                                <factor>
## GSM1012777 Next, we ran RSEM v1.11 with default parameters on these alignments to estimate expression levels. RSEM’s gene level expression estimates (tau) were multiplied by 1,000,000 to obtain transcript per million (TPM) estimates for each gene.
## GSM1012778 Next, we ran RSEM v1.11 with default parameters on these alignments to estimate expression levels. RSEM’s gene level expression estimates (tau) were multiplied by 1,000,000 to obtain transcript per million (TPM) estimates for each gene.
##            data_processing.2
##                     <factor>
## GSM1012777 Genome_build: mm9
## GSM1012778 Genome_build: mm9
##                                                                                                                                                                                                                                                                                                        data_processing.3
##                                                                                                                                                                                                                                                                                                                 <factor>
## GSM1012777 Supplementary_files_format_and_content: File allGenesTPM.txt represents a matrix of gene expression estimates across all non-MolecularBarcode samples. File umbExp.txt represents a matrix of gene expression estimates across all MolecularBarcode samples.  Linked as supplementary files on Series record.
## GSM1012778 Supplementary_files_format_and_content: File allGenesTPM.txt represents a matrix of gene expression estimates across all non-MolecularBarcode samples. File umbExp.txt represents a matrix of gene expression estimates across all MolecularBarcode samples.  Linked as supplementary files on Series record.
##            platform_id  contact_name        contact_email contact_phone
##               <factor>      <factor>             <factor>      <factor>
## GSM1012777    GPL13112 Rahul,,Satija rsatija@nygenome.org    6177022468
## GSM1012778    GPL13112 Rahul,,Satija rsatija@nygenome.org    6177022468
##            contact_laboratory      contact_institute
##                      <factor>               <factor>
## GSM1012777         Satija Lab New York Genome Center
## GSM1012778         Satija Lab New York Genome Center
##                       contact_address  contact_city contact_state
##                              <factor>      <factor>      <factor>
## GSM1012777 101 Avenue of the Americas New York City            NY
## GSM1012778 101 Avenue of the Americas New York City            NY
##            contact_zip.postal_code contact_country data_row_count
##                           <factor>        <factor>       <factor>
## GSM1012777                   10013             USA              0
## GSM1012778                   10013             USA              0
##               instrument_model library_selection library_source
##                       <factor>          <factor>       <factor>
## GSM1012777 Illumina HiSeq 2000              cDNA transcriptomic
## GSM1012778 Illumina HiSeq 2000              cDNA transcriptomic
##            library_strategy
##                    <factor>
## GSM1012777          RNA-Seq
## GSM1012778          RNA-Seq
##                                                       relation
##                                                       <factor>
## GSM1012777 SRA: http://www.ncbi.nlm.nih.gov/sra?term=SRX190719
## GSM1012778 SRA: http://www.ncbi.nlm.nih.gov/sra?term=SRX190720
##                                                               relation.1
##                                                                 <factor>
## GSM1012777 BioSample: http://www.ncbi.nlm.nih.gov/biosample/SAMN01737621
## GSM1012778 BioSample: http://www.ncbi.nlm.nih.gov/biosample/SAMN01737622
##                                                                             supplementary_file_1
##                                                                                         <factor>
## GSM1012777 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190719
## GSM1012778 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX190/SRX190720
```

