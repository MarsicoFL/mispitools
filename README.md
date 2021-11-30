## Installation

The goal of mispitools is to bring a simulation framework for decision
making in missing person identification cases. You can install the
released version of mispitools from [CRAN](https://CRAN.R-project.org)
with:

## Example

This is an example based on a grandchild identification, first you
should do the simulations:

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

    ## Loading required package: pedtools

Once obtained, false postive (FPR) and false negative rates (FNR) could
be computed. This allows to calculate Matthews correlation coefficient
for a specific LR threshold (T):

    ## [1] "FNR = 0.876 ;  FPR = 0.006 ;  MCC = 0.239325804722877"

Likelihoold ratio distributions under both hypothesis, relatedness and
unrelatedness could be plotted.

![](README_files/figure-markdown_github/LRdist-1.png)

Decision plot brings the posibility of analyzing FPR and FNR for each LR
threshold.

![](README_files/figure-markdown_github/deplot-1.png)

Decision threshold could be calculated. For further reading please see
DOI: 10.1016/j.fsigen.2021.102519

    ## [1] "Decision threshold is: 5"

Please cite this tool as: Marsico, F. L., Vigeland, M. D., Egeland, T.,
& Piñero, M. H. (2021). Making decisions in missing person
identification cases with low statistical power. Forensic science
international: genetics, 102519.
