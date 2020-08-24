Vaginal Transfer and development
================

In this repo, the data and code for producing the results of the
manuscript *Stability of the vaginal microbiota during pregnancy and its
importance for early infant colonization* by Mortensen et al (submitted
to eLife - 2020).

The *processed data* is organized in a phyloseq object and can be
downloaded directly from this repo. 

``` r
load('COPSACbirthmicrobiome_ASV.RData')
```

The code for running the analysis is found in the markdown document
*FullAnalysis.md* including costumized functions. Additionally,
three files; getTransferStats.R, getWinnerStats.R and
inferenceTransferStat.R are used for the analysis.

These are also build into an R-package *MBtransfeR*

``` r
devtools::install_github('mortenarendt/MBtransfeR')
library(MBtransfeR)
```

The **Raw sequence data** is to be found on XX XX.

``` r
phy1 <- subset_samples(phyX, Type=='V' & Time == '36' & DELIVERY=='Normal')
phy2 <- subset_samples(phyX, Type=='F' & Time == '1w' & DELIVERY=='Normal')
```

## Calculate individual transferstats

Here we use only a few permutations as this is time consuming. However,
if you want to get proper inference set nperm to at least 1000.

``` r
res <- randpermutationTransferStats(phy1, phy2, 'dyadnb', nperm = 3)
```
