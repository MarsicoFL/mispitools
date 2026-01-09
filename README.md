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

## The Likelihood Ratio Framework

In missing person cases, we evaluate evidence under two competing hypotheses:

- **H1**: The person of interest (POI) **is** the missing person (MP)
- **H2**: The POI is **not** the MP (comes from the reference population)

The likelihood ratio quantifies how the evidence updates our belief:

$$LR = \frac{P(\text{Evidence} \mid H_1)}{P(\text{Evidence} \mid H_2)}$$

An LR > 1 indicates the evidence is more probable under H1; an LR < 1 favors H2. The LR is **not** a probability of identification—it measures the relative support provided by the evidence.

## Tutorial: Simulating LR Distributions

### Step 1: Genetic Evidence

First, we simulate LR distributions from DNA evidence. This requires defining a pedigree structure connecting the MP to a reference individual.

```r
library(mispitools)
library(forrel)
library(pedtools)

# Define pedigree: grandparent-grandchild relationship
ped <- linearPed(2)  # 3-generation pedigree
ped <- setMarkers(ped, locusAttributes = NorwegianFrequencies[1:15])
ped <- profileSim(ped, N = 1, ids = 2)  # Simulate reference profile

# Simulate LRs: under H1 (POI is MP) and H2 (POI is unrelated)
lr_dna <- sim_lr_genetic(ped, missing = 5, numsims = 500)
lr_dna_df <- lr_to_dataframe(lr_dna)

head(lr_dna_df)
#>      Related  Unrelated
#> 1  1247.3201  0.0023415
#> 2   892.1547  0.0001823
#> ...
```

The `Related` column contains LRs simulated under H1, while `Unrelated` contains LRs under H2.

### Step 2: Non-Genetic Evidence

Non-genetic evidence (sex, age, anthropological features) also contributes to identification:

```r
# Simulate LR distributions for sex and age
lr_sex <- sim_lr_prelim("sex", numsims = 500)
lr_age <- sim_lr_prelim("age", numsims = 500)

head(lr_sex)
#>    Related Unrelated
#> 1    1.863    0.1052
#> 2    1.863    1.8627
#> ...
```

### Step 3: Combining Evidence

Under conditional independence, LRs from different evidence sources multiply:

```r
# Combine DNA + sex + age
lr_total <- lr_combine(lr_dna_df, lr_sex)
lr_total <- lr_combine(lr_total, lr_age)

# Visualize the combined LR distribution
plot_lr_distribution(lr_total)
```

<p align="center">
<img src="man/figures/combined_evidence.png" width="550">
</p>

The separation between distributions under H1 (blue) and H2 (red) reflects the discriminating power of the combined evidence. Greater separation means better ability to distinguish between the two hypotheses.

### Step 4: Database Search Application

In practice, we search databases containing unidentified individuals. The LR ranks candidates by evidential support:

<p align="center">
<img src="man/figures/database_search.png" width="550">
</p>

Candidates are ranked by their combined LR. The individual who is actually the MP (blue) rises to the top of the ranking, demonstrating the practical utility of combining multiple evidence sources.

### Step 5: Decision Analysis

To make decisions, we can compute error rates at different LR thresholds:

```r
# Find optimal threshold balancing false positives and false negatives
threshold <- decision_threshold(lr_total, weight = 10)

# Examine error rates at this threshold
threshold_rates(lr_total, threshold)
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `sim_lr_genetic()` | Simulate LRs from DNA evidence |
| `sim_lr_prelim()` | Simulate LRs from non-genetic evidence |
| `lr_combine()` | Combine independent evidence sources |
| `lr_to_dataframe()` | Convert genetic LR results to data frame |
| `decision_threshold()` | Find optimal classification threshold |
| `threshold_rates()` | Compute error rates at a given threshold |
| `plot_lr_distribution()` | Visualize LR distributions |
| `mispitools_app()` | Interactive Shiny application |

## Citation

Marsico FL, Vigeland MD, Herrera Pinero F, Egeland T (2021). "Making decisions in missing person identification cases with low statistical power." *Forensic Science International: Genetics*, 52, 102519. https://doi.org/10.1016/j.fsigen.2021.102519

## Related Packages

- [forrel](https://github.com/magnusdv/forrel): Forensic pedigree analysis
- [pedtools](https://github.com/magnusdv/pedtools): Pedigree manipulation

## Authors

**Franco L. Marsico** — Creator and Head Maintainer
[![GitHub](https://img.shields.io/badge/GitHub-MarsicoFL-blue?logo=github)](https://github.com/MarsicoFL)

**Main contributors:**
- Suisei Nakagawa [![GitHub](https://img.shields.io/badge/GitHub-SuiseiNakagawa-blue?logo=github)](https://github.com/SuiseiNakagawa)
- Undral Ganbaatar [![GitHub](https://img.shields.io/badge/GitHub-undralg-blue?logo=github)](https://github.com/undralg)

## License

GPL-3
