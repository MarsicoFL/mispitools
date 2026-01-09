<img src="inst/www/MispiIcon.png" align="left" width="120">

# mispitools: Likelihood Ratios in Forensic Sciences

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mispitools)](https://CRAN.R-project.org/package=mispitools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/mispitools?color=blue)](https://cran.r-project.org/package=mispitools)
<!-- badges: end -->

<br clear="left"/>

## Overview

**mispitools** computes likelihood ratios (LRs) for forensic identification, combining genetic (DNA) and non-genetic evidence in missing person investigations.

## Installation

```r
# From CRAN
install.packages("mispitools")

# Development version
devtools::install_github("MarsicoFL/mispitools")
```

## The Problem

In missing person cases, investigators must search databases of unidentified individuals. The challenge: **how to rank candidates and quantify the weight of evidence?**

The likelihood ratio provides a principled answer:

$$LR = \frac{P(\text{Evidence} \mid \text{H}_1: \text{Match})}{P(\text{Evidence} \mid \text{H}_2: \text{No match})}$$

## Practical Example: Database Search

Consider searching a database of 100 candidates using combined evidence:
- **DNA**: 15 STR markers (grandparent-grandchild relationship)
- **Non-genetic**: Sex and age

```r
library(mispitools)
library(forrel)

# Genetic evidence
ped <- linearPed(2)
ped <- setMarkers(ped, locusAttributes = NorwegianFrequencies[1:15])
ped <- profileSim(ped, N = 1, ids = 2)
lr_dna <- sim_lr_genetic(ped, missing = 5, numsims = 500)
lr_dna_df <- lr_to_dataframe(lr_dna)

# Non-genetic evidence
lr_sex <- sim_lr_prelim("sex", numsims = 500)
lr_age <- sim_lr_prelim("age", numsims = 500)

# Combine all evidence
lr_total <- lr_combine(lr_dna_df, lr_sex)
lr_total <- lr_combine(lr_total, lr_age)

# Visualize
plot_lr_distribution(lr_total)
```

### LR Distributions: Separating Matches from Non-Matches

<p align="center">
<img src="man/figures/combined_evidence.png" width="750">
</p>

The combined evidence creates separation between the LR distributions under H1 (true matches) and H2 (non-matches). This separation is what enables effective database searching.

### Database Search: Finding the True Match

<p align="center">
<img src="man/figures/database_search.png" width="750">
</p>

When candidates are ranked by their combined LR, the true match rises to the top. This illustrates the practical utility of combining multiple evidence sources.

## Main Functions

| Function | Purpose |
|----------|---------|
| `sim_lr_genetic()` | Simulate LRs from DNA evidence |
| `sim_lr_prelim()` | Simulate LRs from non-genetic evidence |
| `lr_combine()` | Combine independent evidence sources |
| `lr_to_dataframe()` | Convert genetic LR results |
| `decision_threshold()` | Find optimal classification threshold |
| `plot_lr_distribution()` | Visualize LR distributions |
| `mispitools_app()` | Interactive Shiny application |

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
