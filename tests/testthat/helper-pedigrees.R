## Canonical pedigrees used by cross-engine verification tests.
## Each builder returns a `pedtools::ped` with default labels.
## Referenced from test-vs-*.R files in F1.3, F1.7, F2.6, F3.3, F4.5, F5.3.

ped_parent_child <- function() {
  pedtools::nuclearPed(1)
}

ped_full_sibs <- function() {
  pedtools::nuclearPed(2)
}

ped_half_sibs <- function() {
  pedtools::halfSibPed()
}

ped_grandparent_grandchild <- function() {
  pedtools::linearPed(2)
}

ped_avuncular <- function() {
  pedtools::avuncularPed()
}

## Trim a population frequency vector (from `data(Argentina)`) to its
## top `top_k` alleles, renormalize, and order names by numeric value
## so that both pedtools::marker() and marker_model() index alleles
## identically. Trimming bounds the brute-force state space of the
## R-reference engine (F1.2); the C++ engine (F2.x) and pedprobr
## peeling have no such constraint.
top_k_freqs <- function(freq_table, marker_name, top_k = 4L) {
  vec <- freq_table[[marker_name]]
  alleles <- as.character(freq_table$Allele)
  keep <- vec > 0
  f <- vec[keep]
  a <- alleles[keep]
  ord <- order(f, decreasing = TRUE)
  k <- min(top_k, length(f))
  f <- f[ord][seq_len(k)]
  a <- a[ord][seq_len(k)]
  ord2 <- order(suppressWarnings(as.numeric(a)))
  f <- f[ord2] / sum(f)
  names(f) <- a[ord2]
  f
}
