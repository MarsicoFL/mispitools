## Mispitools: Missing Person Identification Tools

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/mispitools)](https://CRAN.R-project.org/package=mispitools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/mispitools?color=blue)](https://cran.r-project.org/package=mispitools)

<!-- badges: end -->


## Installation

The goal of mispitools is to bring a simulation framework for decision
making in missing person identification cases. You can install it from CRAN typing on your R command line the line presented below:

``` r
install.packages("mispitools")
library(mispitools)
```

You can install too the
versions under development (unstable) of mispitools from [Github](https://github.com/MarsicoFL/mispitools/)
with:
``` r
install.packages("devtools")
library(devtools)
install_github("MarsicoFL/mispitools")
library(mispitools)
```

## Example

This is an example based on a grandchild identification, first you
should do the simulations:

``` r
library(mispitools)
library(forrel)
x = linearPed(2)
x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
x = profileSim(x, N = 1, ids = 2)[[1]]
datasim = makeLRsims(x, missing = 5, 1000, 123)
```

Once obtained, false postive (FPR) and false negative rates (FNR) could
be computed. This allows to calculate Matthews correlation coefficient
for a specific LR threshold (T):

``` r
Trates(datasim, 10)
```

    ## [1] "FNR = 0.678 ;  FPR = 0.018 ;  MCC = 0.404650499729402"

Likelihoold ratio distributions under both hypothesis, relatedness and
unrelatedness could be plotted. 

``` r
LRdist(datasim)
```


Decision plot brings the posibility of analyzing FPR and FNR for each LR threshold. 
It could be obtained doing:

``` r
deplot(datasim)
```

![](README_files/figure-markdown_github/deplot-1.png)

Decision threshold could be calculated. For further reading please see
DOI: 10.1016/j.fsigen.2021.102519

``` r
DTsim(datasim, 10)
```

    ## [1] "Decision threshold is: 6"

It is possible work with preliminary investigation data. We use as example the case of the Abuelas de Plaza de Mayo. It consist in the search of abducted children during the last dictatorship in Argentina. Some function are presented below:

``` r
makePOIprelim()

Output:
    POI-ID        DBD Gender       Birth-type       Birth place
1        1 1977-08-25 female   hospital birth      Buenos Aires
2        2 1979-05-09 female   hospital birth      Buenos Aires
3        3 1976-10-10   male   hospital birth North west region
4        4 1978-04-09 female       home birth           Litoral
5        5 1979-12-11   male   hospital birth              Cuyo
6        6 1978-08-04 female   hospital birth      Buenos Aires
7        7 1977-01-31 female   hospital birth      Buenos Aires
8        8 1976-09-30   male   hospital birth      Buenos Aires
9        9 1979-05-22   male   hospital birth      Buenos Aires
10      10 1978-06-29   male       home birth         Patagonia
...
```

It generates a database of preliminary investigation data. Some features such as gender, declared birth date (DBD), region and type of birth are simulated for a set of persons of interest. Different scenarios of searches could be selected. Also, makeMPprelim() simulates preliminary investigation data for MPs. For instance:

``` r
LRprelim()


```
Could be used to compute a likelihood ratio based on birth dates of the missing person (actual birth date or ABD) and for the person of interest (DBD). Methods for LR computations are described here.



Please cite this tool as: Marsico, F. L. et al(2021). Making decisions in missing person
identification cases with low statistical power. Forensic science
international: genetics, 102519.
