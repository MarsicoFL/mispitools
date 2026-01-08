<img src="inst/www/MispiIcon.png" align="left" width="100">

# mispitools: Likelihood ratios in forensic Sciences

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mispitools)](https://CRAN.R-project.org/package=mispitools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/mispitools?color=blue)](https://cran.r-project.org/package=mispitools)
<!-- badges: end -->

## Overview

**mispitools** is an R package for computing likelihood ratios (LRs) in missing person identification cases. It implements a Bayesian framework for evaluating both genetic (DNA) and non-genetic evidence (biological sex, age, hair color, birth date).

## Theoretical Framework

### The Bayesian Paradigm

Evidence evaluation follows the odds form of Bayes' theorem:

```
Posterior Odds = Prior Odds × LR
```

Where:
- **Prior Odds** = P(H1)/P(H2) before observing evidence
- **LR** = P(E|H1)/P(E|H2) = weight of evidence
- **Posterior Odds** = P(H1|E)/P(H2|E) after observing evidence

The hypotheses are:
- **H1** (prosecution): The unidentified person is the missing person
- **H2** (defense): The unidentified person is not the missing person

### Likelihood Ratio Interpretation

The LR quantifies evidential weight, not probability of identity:

| LR | log₁₀(LR) | Interpretation |
|----|-----------|----------------|
| 1-10 | 0-1 | Weak support for H1 |
| 10-100 | 1-2 | Support for H1 |
| 100-10000 | 2-4 | Strong support for H1 |
| >10000 | >4 | Very strong support for H1 |

For LR < 1, the evidence supports H2 (the person is NOT the missing individual). LR = 1 means uninformative evidence.

### Evidence Combination

Under conditional independence given each hypothesis:

```
LR_total = LR₁ × LR₂ × ... × LRₙ
```

**Critical assumption**: This requires that the evidence sources are conditionally independent. For example, sex and age are typically independent given identity, but hair color and skin color are correlated through genetics.

## Installation

From CRAN:
```r
install.packages("mispitools")
```

Development version:
```r
devtools::install_github("MarsicoFL/mispitools")
```

## Worked Examples

### Example 1: Sex Evidence

**Scenario**: Missing person (MP) is female. Unidentified person observed as female.

**Model**:
- P(observe female | H1, MP=F) = 1 - ε = 0.95 (allowing 5% observation error)
- P(observe female | H2) = 0.50 (population frequency)

```r
library(mispitools)

# Direct calculation
eps <- 0.05       # observation error rate
freq_F <- 0.50    # female frequency in population

LR <- (1 - eps) / freq_F
# LR = 0.95 / 0.50 = 1.9
```

**Result**: LR = 1.9

**Interpretation**: The observation of female sex is 1.9 times more probable if the unidentified person is the missing person than if she is a random person from the population. This provides limited support for H1.

### Example 2: Age Evidence

**Scenario**: MP estimated age 40 years (±6 years). Observed age falls within range.

**Model** (uniform reference):
- P(in range | H1) = 1 - ε = 0.95
- P(in range | H2) = range_width / max_age = 12/80 = 0.15

```r
eps <- 0.05
range_width <- 12  # 2 × 6 years
max_age <- 80

p_range_H2 <- range_width / max_age
LR <- (1 - eps) / p_range_H2
# LR = 0.95 / 0.15 = 6.33
```

**Result**: LR = 6.33

**Interpretation**: Age matching provides moderate support for H1.

### Example 3: Region Evidence

**Scenario**: Search across 6 geographic regions. MP and observed person from same region.

**Model**:
- P(match | H1) = 1 - ε = 0.95
- P(match | H2) = 1/6 ≈ 0.167 (uniform)

```r
eps <- 0.05
nreg <- 6

LR <- (1 - eps) / (1 / nreg)
# LR = 0.95 / 0.167 = 5.7
```

**Result**: LR = 5.7

### Example 4: Combining Evidence

**Scenario**: All three evidence types match (sex, age, region).

```r
LR_sex <- 1.9
LR_age <- 6.33
LR_region <- 5.7

LR_combined <- LR_sex * LR_age * LR_region
# LR_combined = 68.6
# log10(LR) = 1.84
```

**Result**: Combined LR = 68.6 (log₁₀ = 1.84)

**Interpretation**: The combined non-genetic evidence provides moderate support for H1. This can be further combined with DNA evidence if available.

### Example 5: Simulation and Decision Analysis

```r
# Simulate LR distributions under both hypotheses
set.seed(123)
lr_sex <- sim_lr_prelim("sex", numsims = 1000)
lr_age <- sim_lr_prelim("age", numsims = 1000)

# Combine
lr_combined <- lr_combine(lr_sex, lr_age)

# Distribution summaries
# Under H1: median log10(LR) ≈ 0.28
# Under H2: median log10(LR) ≈ -1.00

# Find decision threshold
# weight = 10 means false positives are 10× worse than false negatives
threshold <- decision_threshold(lr_combined, weight = 10)
# threshold = 1.9

# Error rates at this threshold
rates <- threshold_rates(lr_combined, threshold)
# FPR = 0.000 (no false positives)
# FNR = 0.053 (5.3% false negatives)

# Visualize
plot_lr_distribution(lr_combined)
```

### Example 6: Sensitivity Analysis

Understanding how conclusions depend on parameter choices:

```r
# How does sex LR vary with error rate?
sens <- lr_sensitivity("sex", param = "eps", range = c(0.01, 0.20))

# At eps = 0.01: LR = 1.98
# At eps = 0.10: LR = 1.80
# At eps = 0.20: LR = 1.60

plot(sens$param_value, sens$LR, type = "l",
     xlab = "Error rate", ylab = "LR",
     main = "Sex LR sensitivity to observation error")
```

## Function Reference

### LR Calculation

| Function | Evidence Type | Key Parameters |
|----------|--------------|----------------|
| `lr_sex()` | Biological sex | `MPs`, `eps`, `Ps` |
| `lr_age()` | Age range | `MPa`, `MPr`, `epa` |
| `lr_hair_color()` | Hair color (5 categories) | `MPc`, `epc`, `Pc` |
| `lr_birthdate()` | Birth date discrepancy | `ABD`, `DBD`, `alpha` |
| `lr_pigmentation()` | Combined pigmentation | `df` from simulation |
| `lr_combine()` | Combine sources | Two LR dataframes |
| `lr_sensitivity()` | Parameter sensitivity | `evidence_type`, `param` |

### Simulation

| Function | Purpose |
|----------|---------|
| `sim_lr_prelim()` | Simulate LR distributions for non-genetic data |
| `sim_lr_genetic()` | Simulate LR distributions for DNA evidence |
| `sim_poi_prelim()` | Generate synthetic POI database |
| `sim_reference_pop()` | Generate population with correlated pigmentation |

### Decision Analysis

| Function | Purpose |
|----------|---------|
| `decision_threshold()` | Compute optimal LR threshold |
| `threshold_rates()` | Compute FPR, FNR, TPR, TNR, MCC |
| `plot_decision_curve()` | Plot error rate trade-off |
| `plot_lr_distribution()` | Plot LR distributions |

### Interactive Application

```r
mispitools_app()
```

Provides interactive interface for LR calculations, CPT visualization, and decision analysis.

## Assumptions and Limitations

1. **Conditional Independence**: Evidence combination by multiplication assumes conditional independence given each hypothesis. Violations (e.g., correlated traits) may inflate or deflate combined LRs.

2. **Error Rates**: Default error rates (typically ε = 0.05) are illustrative. Actual rates should be estimated from validation studies.

3. **Population Frequencies**: Results depend on reference population choice. Frequencies from inappropriate populations may bias LRs.

4. **Observation Model**: The package uses simplified error models. Real observation errors may be more complex (e.g., asymmetric, context-dependent).

5. **Dirichlet Model**: Birth date LRs use Dirichlet-based estimation. The alpha parameters should reflect empirical data from solved cases.

## References

Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence in missing person cases." *Forensic Science International: Genetics*, 66, 102891. https://doi.org/10.1016/j.fsigen.2023.102891

Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making decisions in missing person identification cases with low statistical power." *Forensic Science International: Genetics*, 52, 102519. https://doi.org/10.1016/j.fsigen.2021.102519

Balding DJ, Steele CD (2015). *Weight-of-Evidence for Forensic DNA Profiles*. 2nd ed. Wiley.

Kling D, Tillmar AO, Egeland T (2014). "Familias 3-Extensions and new functionality." *Forensic Science International: Genetics*, 13, 121-127. https://doi.org/10.1016/j.fsigen.2014.07.016

## Related Packages

- [forrel](https://github.com/magnusdv/forrel): Forensic pedigree analysis and LR computation
- [pedtools](https://github.com/magnusdv/pedtools): Pedigree construction and manipulation

## Authors

Franco L. Marsico (maintainer)

Contributors: Suisei Nakagawa, Undral Ganbaatar

## License

GPL-3
