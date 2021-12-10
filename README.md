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

Please cite this tool as: Marsico, F. L. et al(2021). Making decisions in missing person
identification cases with low statistical power. Forensic science
international: genetics, 102519.
