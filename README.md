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

The code for running the analysis is found in the markdown document [FullAnalysis.md](https://github.com/mortenarendt/VagTransfer/blob/master/FullAnalysis.md) including costumized functions. Additionally, three files; getTransferStats.R, getWinnerStats.R and
inferenceTransferStat.R are used for the analysis.

These are also build into an R-package [MBtransfeR](https://github.com/mortenarendt/MBtransfeR)

``` r
devtools::install_github('mortenarendt/MBtransfeR')
library(MBtransfeR)
```

The **Raw sequence data** is to be found on XX XX.

## References

Child **fecal samples** has been published here: 

*Stokholm, Jakob, Martin J. Blaser, Jonathan Thorsen, Morten A. Rasmussen, Johannes Waage, Rebecca K. Vinding, Ann-Marie M. Schoos et al. "Maturation of the gut microbiome and risk of asthma in childhood." Nature communications 9, no. 1 (2018): 1-10.*

Child **airway samples** has been published here: 

*Mortensen, Martin Steen, Asker Daniel Brejnrod, Michael Roggenbuck, Waleed Abu Al-Soud, Christina Balle, Karen Angeliki Krogfelt, Jakob Stokholm et al. "The developing hypopharyngeal microbiota in early life." Microbiome 4, no. 1 (2016): 1-12.*

*Thorsen, Jonathan, Morten A. Rasmussen, Johannes Waage, Martin Mortensen, Asker Brejnrod, Klaus Bønnelykke, Bo L. Chawes et al. "Infant airway microbiota and topical immune perturbations in the origins of childhood asthma." Nature communications 10, no. 1 (2019): 1-8.*

Child **fecal and airway samples** has been published here: 

*Hjelmsø, Mathis H., Shiraz A. Shah, Jonathan Thorsen, Morten Rasmussen, Gisle Vestergaard, Martin S. Mortensen, Asker Brejnrod et al. "Prenatal dietary supplements influence the infant airway microbiota in a randomized factorial clinical trial." Nature communications 11, no. 1 (2020): 1-10.*

Some of the **vaginal samples** has been published here: 

*Rasmussen, M. A., J. Thorsen, M. G. Dominguez-Bello, M. J. Blaser, M. S. Mortensen, A. D. Brejnrod, S. A. Shah et al. "Ecological succession in the vaginal microbiota during pregnancy and birth." The ISME Journal (2020): 1-11.*


## Calculate individual transferstats

### Extract two matched sets 

``` r
phy1 <- subset_samples(phyX, Type=='V' & Time == '36' & DELIVERY=='Normal')
phy2 <- subset_samples(phyX, Type=='F' & Time == '1w' & DELIVERY=='Normal')
```

Here we use only a few permutations as this is time consuming. However,
if you want to get proper inference set nperm to at least 1000.

``` r
res <- randpermutationTransferStats(phy1, phy2, 'dyadnb', nperm = 3)
```
