## F6.7 — Cross-test of the non-genetic evidence pipeline (F6.1-F6.5)
## against Egeland & Marsico (2026), "Using all available evidence to
## solve kinship cases" (Int J Legal Med), the paper that motivates the
## non-genetic framework. Source md:
##   papers/md/Using_all_available_evidence_to_solve_kinship_case.md
##
## Numeric anchor: Table 2 (Example 3.3). It is the only fully
## closed-form numeric example in the paper: a binary comparison
## feature with frequency alpha = 0.8 (beta = 0.2), misclassification
## m = 0.1, SE data x1=1,x2=2 against MP values y1=1,y2=2, seven DVI
## assignments a1..a7 with a7 (none identified) as the reference. The
## paper prints LRSE, LRFDE, combinedLR and posterior to 4 decimals.
##   * LRSE  = {5.0625, 0.0625, 1.1250, 0.1250, 0.5000, 4.5000, 1.0000}
##   * LRFDE = {12.5,   12.5,   5,      5,      5,      5,      1}
##   * combinedLR = LRSE * LRFDE  (paper Eq. 9: LR = LRFDE * LRSE)
##   * posterior  = pi_i L_i / sum pi_j L_j with equal priors (Eq. 10)
##
## The paper's transition matrix (Eq. 3/4) is exactly the categorical
## confusion-matrix observation model of nongenetic_feature (F6.1): an
## identified pair (UP i matched to MP sigma(i)) contributes
## P(x_i | y_sigma(i), V=M) / P(x_i) = M[y, x] / freq[x], and an
## unidentified UP contributes P(x_i)/P(x_i) = 1. The across-UP product
## is the paper's Eq. 7 independence assumption, i.e.
## evidence_combine(mode = "independent"). Tolerance 1e-12 against the
## exact rationals; the printed Table 2 values are matched to their
## 4-decimal precision.
##
## NOT a numeric anchor (verifier rule: never loosen a tolerance to
## pass; disclose instead): Example 3.4 (the only cross-feature
## dependency / Markov example). Its denominators come from the Fig. 3
## joint-frequency heat-map, an image absent from the md, and its
## printed figures are not even mutually consistent (0.952/0.370 =
## 2.573, yet the paper prints LR(dep) = 2.44; 0.952/0.330 = 2.885 vs
## printed 2.77 — the displayed numerator/denominator are rounded
## independently of the reported LRs). Example 3.4 is therefore
## cross-tested only at the *mechanism* level: markov_se must (a)
## reduce to independent when transition rows equal the H2 marginal
## (paper Eq. 4, "Vi and Mj unrelated => P(x=s) = alpha_s"), and (b)
## make the H2 mass follow the joint instead of the product of
## marginals, reproducing the paper's qualitative claim that ignoring
## dependence overestimates the LR for positively-correlated common
## combinations and underestimates it for negatively-correlated rare
## ones. The exact C++ markov_se joint is already pinned at 1e-12 in
## test-evidence-combine-cpp.R (F6.5); here we tie it to the paper.

## --- Table 2 fixture (paper parameters) -----------------------------
em_alpha <- 0.8                    # P(feature = 1)
em_beta  <- 1 - em_alpha           # P(feature = 2)
em_m     <- 0.1                    # misclassification probability
## Transition matrix M (Eq. 3): rows = MP true value t, cols = UP
## observed value s, M[t, s] = P(x = s | y = t, V = M). Row-stochastic.
em_M     <- matrix(c(1 - em_m, em_m,
                     em_m,     1 - em_m), nrow = 2, byrow = TRUE)
em_freq  <- c("1" = em_alpha, "2" = em_beta)   # H2 population marginal

## Realized LR of one UP<->MP pair through the unified non-genetic
## engine: a "pigmentation" categorical feature whose error model is
## the paper's transition matrix and whose H2 is the population
## marginal. observed = MP true value y; the LR at the realized UP
## value x is p_h1[x] / p_h2[x] = M[y, x] / freq[x].
em_pair_lr <- function(y, x) {
  f <- nongenetic_feature("pigmentation", observed = as.character(y),
                           model = list(reference = "marginal"),
                           db_or_freqs = em_freq, error = em_M)
  cpt <- mispitools:::ng_cpt_cpp_wrap(f)
  i <- match(as.character(x), cpt$state)
  cpt$p_h1[i] / cpt$p_h2[i]
}

## Degenerate point-mass per-feature LR distribution (a realized,
## non-random scalar LR), so the across-UP / SE x FDE multiplication
## of the paper (Eq. 7, Eq. 9) is exercised through the actual
## evidence_combine kernel rather than re-implemented here.
em_point <- function(v) list(log10_lr = log10(v), p_h1 = 1, p_h2 = 1)

em_combine_indep <- function(...) {
  r <- mispitools:::cpp_evidence_combine(list(...), mode = "independent")
  10^r$log10_lr
}

## SE data (Example 3.3): x1 = 1, x2 = 2; MP values y1 = 1, y2 = 2.
## Assignment -> realized LRSE through engine + independent combine.
## NULL pair = unidentified UP (LR factor 1, omitted from the product).
em_LRSE <- function(p1, p2) {
  fs <- list()
  if (!is.null(p1)) fs <- c(fs, list(em_point(em_pair_lr(p1[1], p1[2]))))
  if (!is.null(p2)) fs <- c(fs, list(em_point(em_pair_lr(p2[1], p2[2]))))
  if (length(fs) == 0L) return(1)
  if (length(fs) == 1L)
    return(do.call(em_combine_indep, c(fs, list(em_point(1)))))
  do.call(em_combine_indep, fs)
}

## a_k = list(V1 pair, V2 pair); pair = c(MP_true_y, UP_observed_x).
## a7 (both unidentified) is the reference.
em_assign <- list(
  a1 = list(c(1, 1), c(2, 2)),   # V1->M1, V2->M2
  a2 = list(c(2, 1), c(1, 2)),   # V1->M2, V2->M1
  a3 = list(c(1, 1), NULL),      # V1->M1, V2->*
  a4 = list(c(2, 1), NULL),      # V1->M2, V2->*
  a5 = list(NULL,    c(1, 2)),   # V1->*,  V2->M1
  a6 = list(NULL,    c(2, 2)),   # V1->*,  V2->M2
  a7 = list(NULL,    NULL)       # reference
)

em_LRSE_table2  <- c(5.0625, 0.0625, 1.1250, 0.1250, 0.5000, 4.5000, 1.0000)
em_LRFDE_table2 <- c(12.5, 12.5, 5, 5, 5, 5, 1)               # paper, given
em_post_table2  <- c(0.6570, 0.0081, 0.0584, 0.0065, 0.0260, 0.2336, 0.0104)

test_that("Table 2 LRSE column reproduced by the non-genetic engine", {
  ## Exact rationals: 0.9/0.8 = 1.125, 0.1/0.8 = 0.125, 0.9/0.2 = 4.5,
  ## 0.1/0.2 = 0.5, products thereof. Reference a7 = 1 by construction.
  got <- vapply(em_assign,
                function(a) em_LRSE(a[[1]], a[[2]]), numeric(1))
  ## 1e-12: the only error is engine binding round-trip on exact
  ## terminating-decimal arithmetic.
  expect_equal(unname(got), em_LRSE_table2, tolerance = 1e-12)
  ## And to the paper's printed 4-decimal precision.
  expect_equal(round(unname(got), 4), em_LRSE_table2, tolerance = 5e-5)
})

test_that("per-feature LR distribution matches the paper transition model (Eq. 3/4)", {
  ## MP true = 1: the full P(x | y=1, V=M) / P(x) distribution must be
  ## exactly {(1-m)/alpha at x=1, m/beta at x=2} with H1 = M[1, ] and
  ## H2 = population marginal.
  f <- nongenetic_feature("pigmentation", observed = "1",
                           model = list(reference = "marginal"),
                           db_or_freqs = em_freq, error = em_M)
  d <- mispitools:::per_feature_lr_dist_cpp_wrap(f)
  o <- order(d$log10_lr)
  d <- d[o, ]
  expect_equal(d$p_h1, c(em_m, 1 - em_m), tolerance = 1e-12)        # M[1, ]
  expect_equal(d$p_h2, c(em_beta, em_alpha), tolerance = 1e-12)     # freq
  expect_equal(10^d$log10_lr,
               c(em_m / em_beta, (1 - em_m) / em_alpha),
               tolerance = 1e-12)
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(d$p_h2), 1, tolerance = 1e-12)
})

test_that("combinedLR = LRSE * LRFDE via evidence_combine (paper Eq. 9)", {
  lrse <- vapply(em_assign,
                 function(a) em_LRSE(a[[1]], a[[2]]), numeric(1))
  ## Eq. 9: LR(a, a*) = LRFDE * LRSE. LRFDE is genetic input taken from
  ## Table 2 as published (out of the F6 non-genetic verifier scope);
  ## the multiplicative composition is the kernel under test.
  combined <- mapply(function(s, g) em_combine_indep(em_point(s),
                                                     em_point(g)),
                     lrse, em_LRFDE_table2)
  expect_equal(unname(combined),
               em_LRSE_table2 * em_LRFDE_table2, tolerance = 1e-12)
  expect_equal(round(unname(combined), 4),
               c(63.2812, 0.7813, 5.6250, 0.6250, 2.5000, 22.5000, 1.0000),
               tolerance = 5e-5)
})

test_that("Table 2 posterior column reproduced (paper Eq. 10, equal priors)", {
  lrse <- vapply(em_assign,
                 function(a) em_LRSE(a[[1]], a[[2]]), numeric(1))
  combined <- mapply(function(s, g) em_combine_indep(em_point(s),
                                                     em_point(g)),
                     lrse, em_LRFDE_table2)
  ## Eq. 10 with equal priors pi_i = 1/N collapses to normalising the
  ## combined LRs (all expressed relative to a7).
  post <- combined / sum(combined)
  expect_equal(round(unname(post), 4), em_post_table2, tolerance = 5e-5)
})

## --- Example 3.4: cross-feature dependency (mechanism only) ----------

test_that("markov_se reduces to independent when transitions = H2 marginal (Eq. 4)", {
  ## Paper Eq. 4: when Vi and Mj are unrelated, P(x = s) = alpha_s, so
  ## the dependency model must collapse to the independence model.
  fa <- list(log10_lr = log10(c(0.9, 0.1) / c(0.5, 0.5)),
             p_h1 = c(0.9, 0.1), p_h2 = c(0.5, 0.5))
  fb <- list(log10_lr = log10(c(0.8, 0.2) / c(0.4, 0.6)),
             p_h1 = c(0.8, 0.2), p_h2 = c(0.4, 0.6))
  Tmarg <- matrix(c(0.4, 0.6,
                    0.4, 0.6), nrow = 2, byrow = TRUE)  # every row = m2
  ms <- mispitools:::cpp_evidence_combine(list(fa, fb), mode = "markov_se",
                                          transition = list(Tmarg))
  ie <- mispitools:::cpp_evidence_combine(list(fa, fb), mode = "independent")
  expect_equal(ms$log10_lr, ie$log10_lr, tolerance = 1e-12)
  expect_equal(ms$p_h1, ie$p_h1, tolerance = 1e-12)
  expect_equal(ms$p_h2, ie$p_h2, tolerance = 1e-12)
})

test_that("markov_se H2 follows the joint; ignoring dependence over/under-estimates (Example 3.4)", {
  ## Two pigmentation features, 2 states each. H1 numerator is the
  ## same for dep and ind (paper). Marginals m1, m2 and an explicit
  ## positively-correlated joint J (J[1,1] > m1[1]*m2[1]) encoded as a
  ## transition T = J / rowSums so that, under markov_se, the H2 mass
  ## of cell (i,j) is m1[i] * T[i,j] = J[i,j], while independent gives
  ## m1[i] * m2[j].
  m1 <- c(0.6, 0.4)                       # hair-like marginal
  m2 <- c(0.7, 0.3)                       # eye-like marginal
  J  <- matrix(c(0.50, 0.10,              # J[1,1] = .50 > .42 (pos. corr)
                 0.20, 0.20), nrow = 2, byrow = TRUE)
  expect_equal(rowSums(J), m1, tolerance = 1e-12)   # consistent margins
  expect_equal(colSums(J), m2, tolerance = 1e-12)
  Tdep <- J / rowSums(J)                  # row-stochastic

  ## Arbitrary positive H1 (the "match" numerator is shared by dep/ind,
  ## so it cancels in the dep/ind ratio exactly as in the paper).
  h1a <- c(0.95, 0.05); h1b <- c(0.95, 0.05)
  fa <- list(log10_lr = log10(h1a / m1), p_h1 = h1a, p_h2 = m1)
  fb <- list(log10_lr = log10(h1b / m2), p_h1 = h1b, p_h2 = m2)

  dep <- mispitools:::cpp_evidence_combine(list(fa, fb), mode = "markov_se",
                                           transition = list(Tdep))
  ind <- mispitools:::cpp_evidence_combine(list(fa, fb),
                                           mode = "independent")

  ## Total H1 / H2 mass preserved under both models.
  expect_equal(sum(dep$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(dep$p_h2), 1, tolerance = 1e-12)
  expect_equal(sum(ind$p_h1), 1, tolerance = 1e-12)

  ## Realized "common match" cell (1,1): denominator is J[1,1] under
  ## the dependency model and m1[1]*m2[1] under independence.
  num   <- h1a[1] * h1b[1]
  lrdep <- num / J[1, 1]
  lrind <- num / (m1[1] * m2[1])
  ## Positive correlation => independence overestimates the LR
  ## (paper Example 3.4, "ignoring dependence -> overestimation").
  expect_gt(lrind, lrdep)
  expect_equal(lrind / lrdep, J[1, 1] / (m1[1] * m2[1]),
               tolerance = 1e-12)

  ## Negatively-correlated rare cell (1,2): J[1,2] < m1[1]*m2[2]
  ## => independence underestimates (paper's second case).
  expect_lt(J[1, 2], m1[1] * m2[2])
  lrdep2 <- num / J[1, 2]
  lrind2 <- num / (m1[1] * m2[2])
  expect_lt(lrind2, lrdep2)

  ## The engine's H2 mass for the dependency model must literally be
  ## the joint J at the matching support point, and the product of
  ## marginals under independence — the discriminating fact of the
  ## paper's whole "Modelling conditional dependency" section.
  ##   most-extreme log10_lr row = the (1,1) common cell.
  i_dep <- which.max(dep$log10_lr)
  i_ind <- which.max(ind$log10_lr)
  expect_equal(dep$p_h2[i_dep], J[1, 1],            tolerance = 1e-12)
  expect_equal(ind$p_h2[i_ind], m1[1] * m2[1],      tolerance = 1e-12)
  expect_equal(dep$p_h1[i_dep], num,                tolerance = 1e-12)
})
