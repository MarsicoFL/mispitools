<img src="inst/www/MispiIcon.png" align="left" width="120">

# mispitools: Likelihood Ratios in Forensic Sciences

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mispitools)](https://CRAN.R-project.org/package=mispitools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/mispitools?color=blue)](https://cran.r-project.org/package=mispitools)
<!-- badges: end -->

<br clear="left"/>

## Overview

**mispitools** is an R package for computing likelihood ratios (LRs) in forensic identification cases, particularly missing person investigations. It provides a rigorous Bayesian framework for evaluating both genetic (DNA) and non-genetic evidence.

## Installation

```r
# From CRAN
install.packages("mispitools")

# Development version
devtools::install_github("MarsicoFL/mispitools")
```

## Theoretical Framework

### The Bayesian Paradigm

Evidence evaluation in forensic science follows the odds form of Bayes' theorem:

```
Posterior Odds = Prior Odds × LR
```

Where:
- **Prior Odds** = P(H1)/P(H2) before observing evidence
- **LR** = P(E|H1)/P(E|H2) = weight of evidence
- **Posterior Odds** = P(H1|E)/P(H2|E) after observing evidence

The competing hypotheses are:
- **H1** (prosecution): The unidentified person IS the missing person
- **H2** (defense): The unidentified person is NOT the missing person

### The Likelihood Ratio

The LR quantifies the evidential weight—how much more (or less) probable the evidence is under H1 compared to H2:

$$LR = \frac{P(E \mid H_1)}{P(E \mid H_2)}$$

**Interpretation:**
- LR > 1: Evidence supports H1 (identification)
- LR < 1: Evidence supports H2 (exclusion)
- LR = 1: Evidence is uninformative

| log₁₀(LR) | Interpretation |
|-----------|----------------|
| 0 to 1 | Weak support for H1 |
| 1 to 2 | Moderate support for H1 |
| 2 to 4 | Strong support for H1 |
| > 4 | Very strong support for H1 |

### Evidence Combination

When multiple independent evidence sources are available, the combined LR is the product of individual LRs:

$$LR_{total} = LR_1 \times LR_2 \times \cdots \times LR_n$$

**Critical assumption**: This requires conditional independence of evidence given each hypothesis. For example, sex and age are typically independent, but hair color and skin color may be correlated through genetics.

## Quick Example

```r
library(mispitools)

# Simulate LR distributions for sex and age evidence
lr_sex <- sim_lr_prelim("sex", numsims = 1000)
lr_age <- sim_lr_prelim("age", numsims = 1000)

# Combine evidence
lr_combined <- lr_combine(lr_sex, lr_age)

# Visualize the discriminating power
plot_lr_distribution(lr_combined)

# Find optimal decision threshold
decision_threshold(lr_combined, weight = 10)
```

## LR Distributions: H1 vs H2

The following plot illustrates the core concept: when combining sex and age evidence, the LR distributions under H1 (the person IS the missing individual) and H2 (the person is NOT) show clear separation. This separation represents the discriminating power of the evidence.

<p align="center">
<img src="man/figures/lr_distribution.png" width="700">
</p>

Greater separation between distributions indicates stronger evidence. The overlap region represents cases where misclassification is possible, highlighting the importance of combining multiple evidence sources.

## Decision Analysis

Choosing a decision threshold involves balancing two types of errors:
- **False Positives (FP)**: Declaring a match when there is none
- **False Negatives (FN)**: Missing a true match

<p align="center">
<img src="man/figures/decision_curve.png" width="700">
</p>

In forensic contexts, false positives are typically considered more serious than false negatives. The `decision_threshold()` function finds the optimal threshold based on the relative weight assigned to each error type.

## Main Functions

| Function | Purpose |
|----------|---------|
| `sim_lr_prelim()` | Simulate LR distributions for non-genetic evidence |
| `sim_lr_genetic()` | Simulate LR distributions for DNA evidence |
| `lr_combine()` | Combine independent evidence sources |
| `decision_threshold()` | Find optimal classification threshold |
| `plot_lr_distribution()` | Visualize H1 vs H2 distributions |
| `mispitools_app()` | Launch interactive Shiny application |

## Interactive Application

```r
mispitools_app()
```

The Shiny app provides an interactive interface for LR calculations, visualization, and decision analysis.

## Citation

Marsico FL, Vigeland MD, et al (2021). "Making decisions in missing person identification cases with low statistical power." *Forensic Science International: Genetics*, 52, 102519. https://doi.org/10.1016/j.fsigen.2021.102519

## Related Packages

- [forrel](https://github.com/magnusdv/forrel): Forensic pedigree analysis
- [pedtools](https://github.com/magnusdv/pedtools): Pedigree manipulation

## Authors

**Franco L. Marsico** — Head Maintainer
[![GitHub](https://img.shields.io/badge/GitHub-MarsicoFL-blue?logo=github)](https://github.com/MarsicoFL)

**Development Team:**
- Suisei Nakagawa [![GitHub](https://img.shields.io/badge/GitHub-SuiseiNakagawa-blue?logo=github)](https://github.com/SuiseiNakagawa)
- Undral Ganbaatar [![GitHub](https://img.shields.io/badge/GitHub-undralg-blue?logo=github)](https://github.com/undralg)

## License

GPL-3
