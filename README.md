<img src="/README_files/figure-markdown_github/MispiIcon.png" align="left" width="110">


## Mispitools: Missing Person Identification Tools

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/mispitools)](https://CRAN.R-project.org/package=mispitools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/mispitools?color=blue)](https://cran.r-project.org/package=mispitools)

<!-- badges: end -->

## About mispitools
'mispitools' is an open source software package written in R statistical language.
It consist in a set of decision making tools to conduct missing person searches. 
Particularly, it allows computing several features, from non-genetic based LRs to optimal LR threshold for declaring potential matches in DNA-based database search.
mispitools imports forrel, https://doi.org/10.1016/j.fsigen.2020.102376, and pedtools packages, https://doi.org/10.1016/C2020-0-01956-0.
More recently 'mispitools' incorporates preliminary investigation data based LRs. Statistical weight of different traces of evidence such as biological sex, age and hair color are presented. 
For citing mispitools please use the following references: Marsico and Caridi, 2023, http://dx.doi.org/10.2139/ssrn.4331033, and  Marsico, Vigeland et al. 2021, https://doi.org/10.1016/j.fsigen.2021.102519.


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

Now you can analyze the mispitools documentation, it has a description for all functions and parameters.

```r 
?mispitools
```


## Computing non-genetic based LRs

Now you are able to compute conditional probability tables, firstly you
can analyze the different parameters from the documentation.

``` r
?CPT_POP
```

You can compute the phenotype probabilities of three variables, sex, age
and hair color in the population. For simplification, the reference ages
distribution in the population are assumed as uniform (the function
could be easily adapted to incorpore a dataset with the specified
frequencies).

``` r
CPT_POP(
  propS = c(0.5, 0.5),
  MPa = 40,
  MPr = 6,
  propC = c(0.3, 0.2, 0.25, 0.15, 0.1))
```

    ##        [,1]  [,2]    [,3]    [,4]   [,5]
    ## F-T1 0.0225 0.015 0.01875 0.01125 0.0075
    ## F-T0 0.1275 0.085 0.10625 0.06375 0.0425
    ## M-T1 0.0225 0.015 0.01875 0.01125 0.0075
    ## M-T0 0.1275 0.085 0.10625 0.06375 0.0425

The obtained matrix represent the probabilities of the phenotypes in the
reference population. Note that in the followin case, the parameters
remains the same, but changin the MP range change the population
probabilities.

``` r
CPT_POP(
  propS = c(0.5, 0.5),
  MPa = 40,
  MPr = 15,
  propC = c(0.3, 0.2, 0.25, 0.15, 0.1))
```

    ##         [,1]   [,2]     [,3]     [,4]    [,5]
    ## F-T1 0.05625 0.0375 0.046875 0.028125 0.01875
    ## F-T0 0.09375 0.0625 0.078125 0.046875 0.03125
    ## M-T1 0.05625 0.0375 0.046875 0.028125 0.01875
    ## M-T0 0.09375 0.0625 0.078125 0.046875 0.03125

In the same way, MP conditioned probabilities could be computed. Again,
we first see the documentation:

``` r
?CPT_MP
```

Then, we can select a specified MP. One of the parameters is epc, that
comes from the function Cmodel(). Lets see that function:

``` r
?Cmodel()
```

It has two options, uniform, that add the same ep for all combinations
of colors, and custom, that specify a specific value for each pair. Here
we select the custom.

``` r
Cmodel(
  errorModel = "custom",
  ep = 0.01,ep12 = 0.01,ep13 = 0.005,
  ep14 = 0.01,ep15 = 0.003,ep23 = 0.01,
  ep24 = 0.003,ep25 = 0.01,ep34 = 0.003,
  ep35 = 0.003,ep45 = 0.01)
```

    ##             [,1]        [,2]        [,3]        [,4]        [,5]
    ## [1,] 0.972762646 0.009727626 0.004863813 0.009727626 0.002918288
    ## [2,] 0.009680542 0.968054211 0.004840271 0.009680542 0.002904163
    ## [3,] 0.004897160 0.009794319 0.979431929 0.002938296 0.002938296
    ## [4,] 0.009746589 0.002923977 0.002923977 0.974658869 0.009746589
    ## [5,] 0.002923977 0.009746589 0.002923977 0.009746589 0.974658869

Now we can specified the phenotype probabilities conditioned on MP
characteristics.

``` r
CPT_MP(MPs = "F", MPc = 1, 
       eps = 0.05, epa = 0.05, 
       epc = Cmodel())
```

    ##                1            2            3            4            5
    ## F-T1 0.877918288 8.779183e-03 4.389591e-03 8.779183e-03 2.633755e-03
    ## F-T0 0.046206226 4.620623e-04 2.310311e-04 4.620623e-04 1.386187e-04
    ## M-T1 0.046206226 4.620623e-04 2.310311e-04 4.620623e-04 1.386187e-04
    ## M-T0 0.002431907 2.431907e-05 1.215953e-05 2.431907e-05 7.295720e-06

Moreover, LR can be computed as follows:

``` r
MP <- CPT_MP(MPs = "F", MPc = 1, 
       eps = 0.05, epa = 0.05, 
       epc = Cmodel())
POP <- CPT_POP(
  propS = c(0.5, 0.5),
  MPa = 40,
  MPr = 6,
  propC = c(0.3, 0.2, 0.25, 0.15, 0.1))

MP/POP
```

    ##                1            2            3            4           5
    ## F-T1 39.01859058 0.5852788586 0.2341115435 0.7803718115 0.351167315
    ## F-T0  0.36240177 0.0054360266 0.0021744106 0.0072480354 0.003261616
    ## M-T1  2.05361003 0.0308041505 0.0123216602 0.0410722006 0.018482490
    ## M-T0  0.01907378 0.0002861067 0.0001144427 0.0003814755 0.000171664

We can see that the sex-age-color: F-T1-1 and M-T1-1 are the only two LR
values over 1, being the former (perfect match) the largest. All these
information could be summarized in the following plot:

``` r
library(ggplot2)
CondPlot(POP,MP)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)<!-- -->

## Calculating DNA-based decision threshold and error rates

In this example, forrel and pedtools packages provides the scafold for
pedigree definition. Allele frequency database from Argentina is used,
provided by mispitools.

``` r
freq = mispitools::getfreqs(Argentina)[1:5]
x = pedtools::linearPed(2)
x = pedtools::setMarkers(x, locusAttributes = freq)
x = forrel::profileSim(x, N = 1, ids = 2)[[1]]
plot(x, hatched = typedMembers(x))
```

![](MispitoolsDNA_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Mispitools allows simulations LR distributions considering both, H1: UP
is MP and H2: UP is not MP, as true

Once obtained, false postive (FPR) and false negative rates (FNR) could
be computed. This allows to calculate Matthews correlation coefficient
for a specific LR threshold (T):

``` r
Trates(datasim, 10)
```

    ## [1] "FNR = 0.87 ;  FPR = 0 ;  MCC = 0.263664022152322"

Likelihoold ratio distributions under both hypothesis, relatedness and
unrelatedness could be plotted.

``` r
LRdist(datasim, type = 2)
```
<img src="README_files/figure-markdown_github/newplot.png" width="10" height="10">


Or other plotting option:

``` r
LRdist(datasim, type = 1)
```
<img src="README_files/figure-markdown_github/New3.png" width="10" height="10">

Decision plot brings the posibility of analyzing FPR and FNR for each LR
threshold. It could be obtained doing:

``` r
deplot(datasim)
```

<img src="README_files/figure-markdown_github/newplot2.png" width="5" height="5">
