# mispitools 2.0.0

This release adds an exact computational engine alongside the simulation
workflow of the 1.x series. Where 1.x estimated likelihood-ratio
distributions by Monte Carlo, the new layer computes them exactly from the
pedigree and the allele frequencies, so quantities such as the expected
weight of evidence or a tail probability no longer carry simulation error.
The kernel is written in C++ and reached through Rcpp.

All functions from 1.x remain exported and continue to work.

## New: model layer

* `marker_model()` builds a validated marker model from a pedigree, a
  marker identifier, allele frequencies, and an optional mutation or
  linkage specification. It has `print()` and validation methods, so
  input errors surface at construction rather than mid-computation.
* `nongenetic_feature()` gives non-genetic evidence (biological sex,
  pigmentation, age, birthdate) the same treatment as a marker: a
  population distribution, an observed value, and an error rate.
* `get_allele_freqs()` converts the bundled population databases into the
  per-marker format the model layer expects.

## New: exact likelihood-ratio distributions

* `lr_distribution()` returns the full distribution of the profile
  `log10` LR under both hypotheses, not a point estimate. `method =
  "exact"` performs a sparse convolution over the per-marker supports;
  `method = "grid"` projects onto a fixed lattice, which preserves total
  mass and the first moment exactly and keeps the computation bounded
  when the exact support would grow past what memory allows.
* `summary()`, `quantile()` and `plot()` methods for the resulting
  object.
* `per_marker_kl()` and `per_marker_kl_profile()` report the
  Kullback-Leibler divergence between the two hypotheses for each
  marker, which quantifies how much discriminating power each marker
  contributes before any data are observed.

## New: mutation and linkage

* Mutation models: equal, stepwise, and the asymmetric model of Dawid
  (2002).
* Two-marker linked joint distributions computed by Elston-Stewart
  peeling, with recombination handled explicitly rather than assumed
  away.

## New: decision quantities and fragility

* Optimal threshold, false-positive and false-negative rates, ROC curve
  and AUC computed on the exact distribution.
* `concentration_index_positive()`, `leave_one_out()`,
  `calibrate_concentration_cutoff()` and `fragility_report()` measure how
  much of the weight of evidence rests on a single marker, and produce a
  reportable statement for the case file.
* `belief_trajectory()`, `binary_belief_trajectory()` and
  `trajectory_metrics()` follow the posterior as evidence accumulates.

## Verification

The engine is checked against independent implementations rather than
against itself: `pedprobr` and `pedmut` for genotype distributions under
mutation, `Familias` and `forrel` for per-marker likelihood ratios, and
analytic results where they exist. These cross-engine comparisons run as
part of the test suite.

## Performance

The kernel uses OpenMP where the toolchain provides it and falls back to
single-threaded execution otherwise; results are identical either way.
The exact engine enumerates joint genotype states, so its cost grows with
pedigree size and with the number of alleles per marker. For large
pedigrees or full profiles, use `method = "grid"` in `lr_distribution()`.

# mispitools 1.4.0

See the CRAN release notes for the 1.x series.
