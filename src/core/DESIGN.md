# `src/core/` — design specification (mispitools 2.0)

**Status:** F2.1 design. Implementation contract for F2.2 onward.
**Scope:** core C++17 kernel of the mispitools 2.0 evaluation engine.
**Targets driven by this kernel:** R package (via `src/rcpp_bindings.cpp`), WASM
build (via `mispitools_2_loop/wasm/bindings.cpp`, F10+), native CLI (post-2.0).

This document is the contract. F2.2..F5.x implement the modules described here,
in the order set by `STATE.md`. F2.6 / F3.3 / F4.5 / F5.3 verify against the R
reference engine in `R/r_ref_cpt.R` + `R/r_ref_per_marker.R` plus the external
oracles (`pedprobr`, `pedmut`, `fbnet`, `Familias`, `forensIT`).

When SCOUT lifts are realised (`mutation_matrix_stepwise` from fbnet,
`mutation_matrix_proportional` from Familias, etc.), `COPYRIGHTS` at the package
root gets a new entry. See `mispitools_2_loop/scout/` for the per-source notes.

---

## 1. Purpose

`core/` is a pure C++17 library that computes, for a fixed pedigree and a fixed
forensic marker:

1. The joint genotype distribution under H1 (pedigree topology) and under H2
   (POI is unrelated, sampled from HWE).
2. From those two distributions: per-marker LR distribution, per-marker
   bidirectional KL, expected `log10 LR` under each hypothesis.
3. Composition over independent markers (LR distribution convolution).
4. Per-feature analogues for non-genetic evidence (sex, age, region,
   pigmentation, …) under the same algebraic interface.
5. Concentration / fragility indices and decision-theoretic
   FPR/FNR/threshold helpers.

The kernel is *functional* (no global state), *thread-safe* (no static
mutable state on the hot path), and *boundary-agnostic* (no R headers, no
Emscripten headers, no I/O).

The kernel does **not** know about R, pedtools, RcppArmadillo, or the
existence of any external system. All I/O, allocation conversions, and error
reporting happen at the binding layer.

---

## 2. Hard portability rules

The kernel must compile cleanly under three toolchains:

| Toolchain | How `src/core/` is built | Used by |
|---|---|---|
| R/Rcpp (gcc, clang, MSVC via Rtools) | `R CMD INSTALL` via `src/Makevars` | R package |
| Emscripten (`emcc` ≥ 3.x) | `mispitools_2_loop/wasm/emcc_build.sh` | WASM web app |
| Stock gcc/clang (≥ gcc 12, clang 16) | `mispitools_2_loop/scripts/test_core_standalone.sh` | CI smoke + native CLI |

Hard constraints — violations are a CRAN-blocker:

1. **No** `#include <Rcpp.h>`, `<RcppArmadillo.h>`, `<R.h>`, `<Rinternals.h>`,
   `<R_ext/...>` in any `src/core/*.h` or `src/core/*.cpp`.
2. **No** `#include <emscripten/...>` or `EM_JS` in `src/core/`.
3. **No** `Rcpp::stop`, `Rcpp::Rcout`, `Rprintf`, `Rcpp::NumericVector`, or
   any `Rcpp::*` symbol inside `src/core/`.
4. **No** `iostream` in production code paths. Debug-only logging goes through
   a single macro `MISPI_LOG(msg)` that expands to nothing in release builds
   (see §6.3).
5. **No** `throw` of `std::exception` subclasses that may cross the binding
   boundary. Errors are returned as values (§5). The binding layer is the
   only place that may convert an error value into `Rcpp::stop` or a JS
   exception.
6. **C++17 only.** No C++20 (`std::span`, concepts, ranges), no C++14-only
   features either except where it doesn't matter. No compiler-specific
   extensions (`__attribute__`, `__forceinline`).
7. **Header guards** `#ifndef MISPITOOLS_CORE_<MODULE>_H ... #endif`. No
   `#pragma once`.
8. **No external dependencies** outside the C++17 STL. RcppArmadillo /
   Armadillo is allowed *only inside `src/rcpp_bindings.cpp` and (later)
   `src/rcpp/*`*. The core uses `std::vector`, `std::array`, plain arrays,
   `std::optional`, `std::variant` where needed.
9. **No OpenMP `#pragma`** in `src/core/` before F7.1. OpenMP integration is
   confined to the binding layer (or, in F7.1, to a clearly delimited
   parallel section under a feature-test macro `MISPI_HAS_OPENMP`).
10. **No raw `new` / `delete`.** RAII only (`std::vector`, `std::unique_ptr`).

Recommended portability practices (soft):

- `static_assert(sizeof(double) == 8, ...)` at one point of entry — IEEE-754
  is required for KL / LR log-space arithmetic.
- Avoid `<random>` engines outside testing. The kernel itself does not
  sample; sampling lives in the binding layer (R `set.seed`, JS `Math.random`,
  or fixed-seed C++ for the CLI).

---

## 3. Naming & style

- `snake_case` for free functions, variables, struct fields.
- `PascalCase` for type names (`Pedigree`, `Marker`, `JointTable`).
- `UPPER_SNAKE_CASE` for compile-time constants (`MISPI_LOG_FLOOR`).
- Member functions on a struct are *not* required to be members — most
  operations are free functions taking the struct by `const&`. Use member
  functions only for trivial accessors (`size`, `empty`).
- Doxygen-style triple-slash `///` comments on every public function in
  every header (`@brief`, `@param`, `@return`, `@complexity`).
- No comments inside `.cpp` files unless explaining a non-obvious algorithmic
  choice (footnote-style, one line max).
- Public API of `core/` = whatever is declared in `core/*.h`. Internal
  helpers live as static-free functions in the `.cpp` or in an
  `mispitools::core::detail` sub-namespace, never declared in headers.

All public symbols live under `namespace mispitools { namespace core { ... } }`.
Bindings alias with `namespace mc = mispitools::core;`.

---

## 4. Module map

The kernel is split into the following modules. Files already exist as
placeholders (F0.4). Each row links a module to the milestone that fleshes
it out and to the R-side reference engine it is validated against.

| Module | Files (`src/core/`) | First content milestone | R-side oracle | Notes |
|---|---|---|---|---|
| `pedigree` | `pedigree.h` `.cpp` | F2.2 | `R/r_ref_cpt.R::build_joint` | Pedigree POD + topological helpers (founders, ancestral order, descendant flag) |
| `marker` | `marker.h` `.cpp` | F2.2 | `R/marker_model.R::validate_freqs` | Marker POD: alleles (indices), `freqs`, optional `repeat_unit_kind`, optional `numeric_labels` for stepwise |
| `mutation_models` | `mutation_models.h` `.cpp` | F2.3 / F2.4 / F5.1 | `R/r_ref_cpt.R::mutation_matrix_R`, `pedmut::mutationMatrix`, `Familias`, `fbnet::getLocusCPT` | `mutation_matrix_equal` (F2.3), `mutation_matrix_stepwise` (F2.3 lift fbnet/pedmut), `mutation_matrix_proportional` (F2.3 lift Familias), `mutation_matrix_asymmetric` (F5.1) |
| `cpt_engine` | `cpt_engine.h` `.cpp` | F2.2 (none), F2.4 (with mutation) | `R/r_ref_cpt.R::cpt_marker_joint_R` | Founder enumeration + Mendelian peeling. Outputs `JointTable` |
| `linkage` | `linkage.h` `.cpp` | F5.2 | `fbnet::buildBN` + 2-marker `getLocusCPT` | `linked_pair_joint` for two markers with `theta` |
| `kl_engine` | `kl_engine.h` `.cpp` | F3.1 | `R/r_ref_per_marker.R::per_marker_kl_R`, `forensIT::perMarkerKLs` (after upstream fix) | Bidirectional KL on shared support; handles `0 log 0 = 0`; floor at `MISPI_LOG_FLOOR` |
| `lr_dist` | `lr_dist.h` `.cpp` | F4.1 / F4.2 / F4.3 | `R/r_ref_per_marker.R::per_marker_lr_dist_R`, `Familias::FamiliasPosterior$LRperMarker`, `forrel::profileSim` Monte Carlo | Per-marker LR distribution; `lr_dist_compose` exact convolution; quantiles/E/var/ROC |
| `nongenetic_lr` | `nongenetic_lr.h` `.cpp` | F6.3 / F6.4 | `R/lr_sex.R`, `R/lr_age.R`, `R/lr_hair_color.R`, `R/lr_pigmentation.R`, `R/lr_birthdate.R` | CPT under H1/H2 for categorical / continuous / date features |
| `evidence_combine` | `evidence_combine.h` `.cpp` | F6.5 | Egeland-Marsico 2026 reproduction script | `independent` (sum of log10 LRs) and `markov_se` (chain) modes |
| `concentration` | `concentration.h` `.cpp` | F8.2 | `R/fragility_tools.R` | Generalised concentration index `C_W^+`, leave-one-feature-out |
| `decision` | `decision.h` `.cpp` | F4.3 | `R/lr_sensitivity.R` and existing decision functions | FPR / FNR / threshold helpers over an LR distribution |

The file `RcppExports.cpp` is *auto-generated* by `Rcpp::compileAttributes()`
and is not part of the kernel — never hand-edit. The binding wrappers live in
`src/rcpp_bindings.cpp` (top-level under `src/`, intentional: see comment in
that file — `compileAttributes` doesn't recurse into subdirs).

---

## 5. Error handling

Errors that originate inside `core/` must never throw across the kernel
boundary. The chosen idiom is C++17-friendly and avoids the
`tl::expected` dependency.

### 5.1 Result type

For functions that compute a value and may fail:

```cpp
// in src/core/result.h (introduced in F2.2 alongside the first real fn)
namespace mispitools { namespace core {

template <typename T>
struct Result {
    std::optional<T> value;
    std::string error;   // empty iff value.has_value()

    bool ok() const noexcept { return value.has_value(); }
    const T& operator*() const { return *value; }
    T& operator*()             { return *value; }
};

template <typename T>
inline Result<T> ok_result(T&& v) {
    return Result<T>{ std::optional<T>{std::move(v)}, std::string{} };
}

inline auto err_result(std::string msg) {
    // helper templated at call site:  return err_result<T>(...);
    return [msg = std::move(msg)](auto tag) {
        return Result<decltype(tag)>{ std::nullopt, msg };
    };
}

} } // namespace
```

Callers do:

```cpp
auto r = cpt_marker_joint(ped, marker, mut);
if (!r.ok()) {
    return Result<JointTable>{ std::nullopt, r.error };  // bubble up
}
const JointTable& jt = *r;
// ...
```

At the binding boundary (`src/rcpp_bindings.cpp`):

```cpp
auto r = mc::cpt_marker_joint(ped, marker, mut);
if (!r.ok()) Rcpp::stop(r.error);
return as_DataFrame(*r);
```

In WASM bindings the error becomes a `throw new Error(error)` JS-side.

### 5.2 Preconditions

Pre-validated input (`marker_model` already validated by `R/marker_model.R`,
WASM bindings revalidate independently) means most defensive checks happen
at the boundary. Inside `core/`, only invariants the kernel itself relies on
get checked, e.g.:

- `mutation_matrix_*` checks `rate ∈ [0, 1)`, `range ∈ (0, 1)`, etc., and
  returns a `Result` with a precise error message — these are user-input
  checks that the binding layer cannot do without re-implementing them.
- `cpt_marker_joint` checks that `freqs.size() == alleles.size()` and that
  `freqs` sums to 1 within `kSumTol` — this is an invariant guard, not a
  user-input check.

A failed invariant ≠ `assert()`. Use `Result` consistently. `assert()` is
reserved for debug-only structural sanity (e.g. "we got here, the table
must be sorted"). Release builds compile out `assert`.

### 5.3 No exceptions

`core/` does not `throw`. STL functions that can throw (`std::vector::at`,
`std::stoi`, `std::regex_*`) are not used on hot paths. Allocation failure
from `std::bad_alloc` is treated as fatal and is not caught — same policy as
the R session (an OOM is a process-level event).

---

## 6. Numerical conventions

### 6.1 Probability sums

Joint tables coming out of `cpt_marker_joint` must satisfy

```
| sum_g P(g) − 1 |  ≤  kSumTol = 1e-12
```

A larger discrepancy is an invariant failure and returns `Result::error`.
Conversion to log-space inside `kl_engine` / `lr_dist` floors zero / very
small probabilities at:

```
constexpr double MISPI_LOG_FLOOR = 1e-300;     // close to DBL_MIN
inline double safe_log(double p) {
    return std::log(p < MISPI_LOG_FLOOR ? MISPI_LOG_FLOOR : p);
}
```

The choice `1e-300` (not `DBL_MIN ≈ 2.22e-308`) leaves a 10⁸ factor of safety
under subsequent additions. Conversion log₁₀ uses `std::log10`.

### 6.2 KL with zero entries

Convention: `0 · log(0 / x) = 0` (limit). Implementation: any term with
`P[g] == 0` contributes 0, *regardless of the value of `Q[g]`*. A term with
`Q[g] == 0` and `P[g] > 0` is an invariant violation under H2 (incompatible
support) and returns `Result::error` from `per_marker_kl`. This case can
only arise with `mutation = none` plus an evidentially incompatible
configuration; documented in the F3 error path.

### 6.3 Logging

The kernel does not write to stdout, stderr, R console, or any sink. The
following macro is the only logging surface:

```cpp
#ifdef MISPI_DEBUG
#  define MISPI_LOG(stream) do { std::cerr << "[core] " << stream << "\n"; } while (0)
#else
#  define MISPI_LOG(stream) do { } while (0)
#endif
```

`MISPI_DEBUG` is set only by the standalone test harness in
`mispitools_2_loop/scripts/test_core_standalone.sh`. R / WASM builds leave it
undefined.

### 6.4 Determinism

For fixed input the kernel produces bit-identical output across runs and
across the three toolchains. No floating-point summation reordering on hot
paths. When OpenMP arrives in F7.1, parallel reductions use Kahan-compensated
summation only if a regression is detected; otherwise the canonical order
(per-marker outer loop) preserves determinism.

---

## 7. Types

### 7.1 `Pedigree`

```cpp
namespace mispitools { namespace core {

using MemberIndex = std::int32_t;     // 0-based, contiguous over the pedigree
constexpr MemberIndex kNoParent = -1;

/// @brief Plain-old-data pedigree representation.
///
/// All vectors have length `n_members`. Indices into other vectors use the
/// dense `MemberIndex` numbering. External labels (the R `id` column) are
/// erased at the binding boundary; core operates on indices only.
///
/// Ordering invariant: members are listed in a topologically-valid order,
/// i.e. for any i, `father[i] < i` and `mother[i] < i` when they are not
/// `kNoParent`. The binding layer enforces this; the kernel may rely on it.
struct Pedigree {
    MemberIndex n_members = 0;
    std::vector<MemberIndex> father;     // size n_members; kNoParent for founders
    std::vector<MemberIndex> mother;     // size n_members; kNoParent for founders
    std::vector<std::uint8_t> sex;       // 0 = unknown, 1 = male, 2 = female
    std::vector<std::uint8_t> is_typed;  // 0/1 per member (this marker)
    std::vector<MemberIndex> founders;   // cached: {i : father[i]==kNoParent && mother[i]==kNoParent}
    std::vector<MemberIndex> nonfounders;// cached: complement of founders, in topo order
    MemberIndex poi = kNoParent;         // person of interest (excluded under H2)
};

}} // namespace
```

Helpers in `pedigree.h`:

```cpp
bool is_founder(const Pedigree& p, MemberIndex i) noexcept;
bool has_descendant_typed(const Pedigree& p, MemberIndex i) noexcept;
std::vector<MemberIndex> ancestral_order_excluding(
    const Pedigree& p, MemberIndex excluded);
```

Note on linkage: the same `Pedigree` is reused across markers in a case.
Linkage data attaches to `Marker` (§7.2), not to `Pedigree`.

### 7.2 `Marker`

```cpp
namespace mispitools { namespace core {

using AlleleIndex = std::int32_t;    // 0..n_alleles-1

/// @brief Single forensic marker with population frequencies and typings.
///
/// `typing[i]` is the observed genotype of member `i` if `is_typed[i] == 1`,
/// encoded as a sorted pair `{a1, a2}` with `a1 <= a2`. When the member is
/// not typed, the typing entry is `{kNoAllele, kNoAllele}`.
struct Genotype {
    AlleleIndex a1 = -1;   // sorted: a1 <= a2; -1 means missing
    AlleleIndex a2 = -1;
};

struct Marker {
    std::string id;                       // e.g. "D3S1358"
    AlleleIndex n_alleles = 0;
    std::vector<double> freqs;            // size n_alleles, sums to 1 within kSumTol
    std::vector<double> numeric_labels;   // size n_alleles, used by stepwise; NaN if not parseable
    std::vector<Genotype> typing;         // size n_members of the bound Pedigree; -1 entries for untyped
    bool always_lumpable = false;         // precomputed (F2.5+) for shortcut paths
};

}} // namespace
```

Allele labels (the string names "12", "13.2", ...) are erased at the
boundary. The kernel works on integer `AlleleIndex` and on the parallel
numeric vector for stepwise step distances.

### 7.3 `MutationModel`

```cpp
namespace mispitools { namespace core {

enum class MutationKind : std::uint8_t {
    None       = 0,   // identity
    Equal      = 1,   // M[i,i] = 1-R, M[i,j]=R/(K-1) for i!=j
    Stepwise   = 2,   // M[i,j] = (R / sum r^|s_i-s_k|) * r^|s_i-s_j|
    Proportional = 3, // M[i,j] = alpha*p[j] off-diag, alpha = R / sum p(1-p)
    Asymmetric = 4    // F5.1; biased stepwise with bias u
};

struct MutationModel {
    MutationKind kind = MutationKind::None;
    double rate  = 0.0;     // R, in [0, 1)
    double range = 0.0;     // r (stepwise/asymmetric), in (0, 1)
    double rate2 = 0.0;     // out-of-microgroup rate, default 0
    double bias  = 0.5;     // u for asymmetric, in [0, 1]
};

}} // namespace
```

Note: paternal and maternal mutation matrices are *separate* in `Marker`
in F5+ (sex-specific mutation). The `MutationModel` above is shared in the
default case; the F5.1 spec introduces an optional `paternal/maternal` pair.

### 7.4 `JointTable` (sparse joint genotype distribution)

The CPT joint over the typed-and-relevant subset of members is sparse — most
combinations of `n_alleles^{2 * n_members}` have zero probability under
Mendelian propagation. Storage:

```cpp
namespace mispitools { namespace core {

using GenotypeIndex = std::int32_t;   // 0..G-1 where G = n_alleles*(n_alleles+1)/2

/// @brief Sparse joint over (genotype-of-each-relevant-member, prob).
///
/// `member_ids[k]` is the dense index inside Pedigree::n_members.
/// Each row in `rows` is a tuple of GenotypeIndex of length n_members_relevant,
/// stored as a flat int32 vector indexed by `row * n_members_relevant + k`.
/// `p_h1[row]` / `p_h2[row]` carry the matched H1 / H2 probabilities.
///
/// Rows are sorted lex on (genotype indices), allowing two JointTables on
/// the same support to be merged in O(rows) by a parallel scan.
struct JointTable {
    std::int32_t n_members_relevant = 0;
    std::vector<MemberIndex> member_ids;            // length n_members_relevant
    std::vector<GenotypeIndex> rows_flat;           // length n_rows * n_members_relevant
    std::vector<double> p_h1;                       // length n_rows
    std::vector<double> p_h2;                       // length n_rows

    std::size_t n_rows() const noexcept { return p_h1.size(); }
    bool empty() const noexcept { return p_h1.empty(); }
};

}} // namespace
```

Helpers (in `cpt_engine.h`):

```cpp
GenotypeIndex pair_to_genotype_index(AlleleIndex a1, AlleleIndex a2, AlleleIndex n_alleles) noexcept;
std::pair<AlleleIndex, AlleleIndex> genotype_index_to_pair(GenotypeIndex g, AlleleIndex n_alleles) noexcept;
inline GenotypeIndex n_genotypes(AlleleIndex n_alleles) noexcept {
    return n_alleles * (n_alleles + 1) / 2;
}
```

Allele pair to genotype index uses the standard upper-triangular numbering:

```
g(a1, a2) = a1 * n_alleles - a1 * (a1 - 1) / 2 + (a2 - a1)     with a1 <= a2
```

The Pedigree-side R reference uses the column-major lex order
`a1 + (a2-1) * G_a2` (1-indexed); the C++ side uses the row-major upper-tri
encoding above. The binding layer maps between them — neither encoding leaks
across the boundary.

### 7.5 `LrDist`

The per-marker LR distribution (`per_marker_lr_dist_R` output):

```cpp
namespace mispitools { namespace core {

/// @brief Sparse LR distribution: (log10 LR, P(.|H1), P(.|H2)) triples.
struct LrDist {
    std::vector<double> log10_lr;   // sorted ascending, aggregated (no duplicates)
    std::vector<double> p_h1;       // same length
    std::vector<double> p_h2;       // same length
    bool has_pos_inf = false;       // a non-empty mass at +Inf (mutation=none, H2 incompatible)
    bool has_neg_inf = false;       // mass at -Inf (rare; only when both supports differ)
};

}} // namespace
```

Convolution (`lr_dist_compose`) takes a `std::vector<LrDist>` and returns a
`LrDist` on the composed grid. Exact convolution is feasible when the total
support of each marker is small (≤ ~50 entries) and number of markers small
(≤ ~30). Beyond that, F4.2 introduces a grid heuristic with controlled error.

---

## 8. Algorithms

### 8.1 `cpt_marker_joint` — F2.2, F2.4

Mirrors `R/r_ref_cpt.R::cpt_marker_joint_R`. The R-reference engine is the
oracle; the C++ engine reproduces its output bit-for-bit (1e-12 tol) and
then beats it on speed.

Algorithm (Elston-Stewart peeling with founder enumeration over informative
founders only):

```
Inputs:  Pedigree P, Marker M, MutationModel mu, MemberIndex poi
Outputs: JointTable joint (with p_h1, p_h2 aligned on the same support)

1. Pre-process:
   1a. Compute G = n_genotypes(M.n_alleles).
   1b. Build HWE prior over genotypes: hwe[g] = p_i^2 (i==j) or 2 p_i p_j (i!=j).
   1c. Build mutation matrix M_mut (K×K) from mu.
   1d. Build transmission table T_trans[g][a] = P(child receives allele a |
       parent has genotype g), applying M_mut as right-multiplication:
       T_trans = T_meiotic %*% M_mut.
   1e. Build child_dist[g_pat][g_mat][g_child] = sum over inheritance
       events under T_trans (vectorised once, reused per child).
   1f. Build child_dist_one_missing[g_known][g_child] =
       sum_{g_other} hwe[g_other] * child_dist[g_known][g_other][g_child].

2. Identify informative founders:
   - founder f is informative iff EXISTS a typed descendant of f.
   - non-informative founders are integrated analytically (they contribute
     a factor of 1; their genotype gets imputed against HWE when descendants
     need it via child_dist_one_missing).

3. Enumerate H1 joint:
   - Start with one row, prob = 1, no member set.
   - For each informative founder f in some fixed order, expand the row set
     by `G` copies, set states[f][r] in {0..G-1}, multiply prob by hwe[g].
   - For each non-founder in topological order:
     - Look up father and mother indices.
     - If both parent genotypes are already in the state vector: multiply
       prob by child_dist[g_pat][g_mat][g_child].
     - If exactly one parent is known: multiply by child_dist_one_missing.
     - Drop rows with prob <= 0 immediately (sparsity gain).
   - Result: JointTable.rows_flat carries (one column per relevant member),
     joint.p_h1[row] = product so far.

4. Enumerate H2 joint:
   - Same as H1 but with POI excluded. POI's genotype is sampled
     independently from HWE; on every typed member that has POI as a parent,
     the parent is treated as "missing" (imputed under HWE).
   - Each H2 row r has joint.p_h2[r] = sum_{g_POI} hwe[g_POI] * (peeled
     sub-pedigree prob), which factorises into hwe[g_POI] times the
     sub-pedigree marginal.

5. Align supports:
   - Expand H2 rows over G choices of POI's genotype to make them
     comparable to H1.
   - Merge into a combined sparse table: row = full state of all members
     (POI included). For rows present in H1 but missing in H2, p_h2 = 0
     (only if H2 truly cannot reach that POI under HWE — should not happen
     for nonzero mutation); for rows in H2 but missing in H1, p_h1 = 0.
   - Drop rows where both p_h1 == 0 and p_h2 == 0.

6. Sort rows lex on (member_0, member_1, ...) to enable later O(rows) joins
   between H1 and H2 in kl_engine / lr_dist.
```

Complexity (worst case): O(G^{F_eff} · n_nonfounders · G²) where `F_eff` is
the number of informative founders. For MP scope (≤ 6 informative founders,
≤ ~15 alleles after lumping, ~5 non-founders), this is ~10⁶..10⁸ ops, well
under the 50 ms wall-clock budget for a single marker.

#### 8.1.1 Reduction techniques

- **Drop-zero rows** after each propagation step (matches R-ref).
- **Allele lumping** (F2.5): non-observed alleles fold into a single "other"
  state when the marker is `always_lumpable` or when an explicit lump
  partition is provided. Exact under Kemeny-Snell lumpability.
- **Sparse transmission tables**: when `n_alleles > 10` and stepwise
  mutation is the model, `T_trans` has structured sparsity — only step-near
  neighbours are non-zero. Stored as a column-CSR for the hot inner loop.
  (Optimisation, F2.5 / F3.4.)

### 8.2 Mutation matrices — F2.3 / F5.1

Direct C++ ports of the R-reference implementations. Lifts identified by
F1.8 SCOUT:

- `mutation_matrix_equal` — SCOUT `pedmut` Lift 1 (3 LOC). Trivial.
- `mutation_matrix_stepwise` — SCOUT `fbnet` Lift 1 (~30 LOC). Already
  validated to 1e-12 against `pedmut::mutationMatrix` and reshaped
  `fbnet::getLocusCPT` in F1.5. Cite fbnet + pedmut in COPYRIGHTS.
- `mutation_matrix_proportional` — SCOUT `Familias` Lift 2 (~10 LOC).
  Reversible w.r.t. `freqs` by construction. Cite Familias.
- `mutation_matrix_asymmetric` — F5.1. Spec to be confirmed under F4.6
  SCOUT pass on Familias bias `u` semantics.

All builders return `Result<std::vector<std::vector<double>>>` so the rare
"undefined model" (negative diagonal under Dawid, `R+R2 > 1` under stepwise,
ill-defined "all stepwise weights zero") propagates cleanly.

Validation helper:

```cpp
/// @brief Check that M is square, finite, entry-wise in [-tol, 1+tol], rowsums == 1 within tol.
/// @complexity O(K^2).
Result<bool> validate_mutation_matrix(
    const std::vector<std::vector<double>>& M,
    double tol = std::sqrt(std::numeric_limits<double>::epsilon()));
```

Note the tightened tolerance vs `pedmut::validateMutationMatrix` (3-decimal,
laxa) — see SCOUT_pedmut.md §"Bugs / caveats".

### 8.3 `per_marker_kl` — F3.1

Mirrors `R/r_ref_per_marker.R::per_marker_kl_R`.

Inputs: a single `JointTable` from `cpt_marker_joint`. Outputs (per marker):

```cpp
struct PerMarkerKL {
    double e_log10_lr_h1 = 0.0;   // E[ log10 LR | H1 ]
    double e_log10_lr_h2 = 0.0;   // E[ log10 LR | H2 ]
    double kl_h1_to_h2  = 0.0;    // KL(P_H1 || P_H2) in nats
    double kl_h2_to_h1  = 0.0;    // KL(P_H2 || P_H1) in nats
    bool   support_mismatch = false;  // H2 zero where H1 positive (incompatible)
};
```

Algorithm: a single linear pass over the sorted rows of the joint, with the
zero-handling convention of §6.2. KL in nats; expectations of log10 LR are
`KL(P||Q)/ln(10)` (with the appropriate sign), avoiding double conversion.

### 8.4 `per_marker_lr_dist` — F4.1

Mirrors `R/r_ref_per_marker.R::per_marker_lr_dist_R`. Pseudocode:

```
Inputs: JointTable joint, optional bool aggregate (default true)
Outputs: LrDist d

for each row r:
  let lr_r = joint.p_h1[r] / joint.p_h2[r]
  if joint.p_h2[r] == 0 && joint.p_h1[r] > 0:
    push (+Inf, p_h1[r], 0) ; d.has_pos_inf = true
  elif joint.p_h1[r] == 0 && joint.p_h2[r] > 0:
    push (-Inf, 0, p_h2[r]) ; d.has_neg_inf = true
  elif joint.p_h1[r] == 0 && joint.p_h2[r] == 0:
    skip
  else:
    push (log10(lr_r), p_h1[r], p_h2[r])

if aggregate:
  sort by log10_lr, then collapse equal-log10_lr entries (sum probs).
```

### 8.5 `lr_dist_compose` — F4.2

Convolution of independent per-marker LR distributions:

```
Inputs: vector<LrDist> per_marker
Output: LrDist d_total

Initialize d_total = { log10_lr: [0.0], p_h1: [1.0], p_h2: [1.0] }
for each pm in per_marker:
  d_new = empty
  for each (lr_old, ph1_old, ph2_old) in d_total:
    for each (lr_pm, ph1_pm, ph2_pm) in pm:
      push (lr_old + lr_pm, ph1_old * ph1_pm, ph2_old * ph2_pm)
  d_total = aggregate(d_new)
```

For ≤ 30 markers with sparse support this is exact and runs under the
performance target (§13). Beyond that, F4.2 falls back to an FFT-style
distribution arithmetic on a fixed log-LR grid (controlled discretisation
error).

### 8.6 `per_feature_*` — F6

Symmetric API for non-genetic features:

```cpp
/// @brief Build (P_H1, P_H2) CPT for a non-genetic feature.
/// @complexity O(card(feature)).
Result<JointTable>
    nongenetic_cpt(const NongeneticFeature& f, const Pedigree& p);

PerMarkerKL per_feature_kl_nongenetic(const JointTable& jt);
LrDist      per_feature_lr_dist_nongenetic(const JointTable& jt);
```

The `NongeneticFeature` struct (defined in F6.1 by `r_api`) plays the role
of `Marker` for the non-genetic side. The same `JointTable` representation
is reused — composition (F6.5) is then identical to §8.5 with both
genetic and non-genetic distributions in the input vector.

### 8.7 `evidence_combine` — F6.5

Two modes:

- `independent`: identical to §8.5.
- `markov_se`: Egeland-Marsico 2026 chain. Spec lives in the paper draft;
  exact formula confirmed when F6.5 is implemented. The kernel exposes both
  through a single function with an enum:

```cpp
enum class CombineMode : std::uint8_t {
    Independent = 0,
    MarkovSE    = 1
};

Result<LrDist> evidence_combine(
    const std::vector<LrDist>& per_feature,
    CombineMode mode,
    const MarkovSEParams& params);   // ignored in Independent mode
```

### 8.8 `concentration` and `decision` — F8.2 / F4.3

Free functions over `LrDist` and over the per-feature vector:

```cpp
double concentration_index_positive(const std::vector<double>& signed_logLRs);
std::vector<double> leave_one_out(const std::vector<double>& signed_logLRs);

struct DecisionResult {
    double threshold;
    double fpr;
    double fnr;
};
DecisionResult choose_threshold(const LrDist& d, double target_fpr);
```

Naming and exact contract finalised when the milestones land.

---

## 9. Sparse representation choice

Alternatives considered and rejected:

| Representation | Pros | Cons | Verdict |
|---|---|---|---|
| Dense N-dim array indexed by all genotypes | Fast indexing, vectorisable | Memory blows up for ≥ 6 members × ≥ 10 alleles (`G^6 = 75⁶ ≈ 10¹¹`) | rejected |
| `std::unordered_map<key, double>` keyed on hashed genotype tuple | Compact, O(1) lookup | Bad cache behaviour, slow iteration, two passes for KL | rejected for hot paths |
| Sorted-on-key `std::vector<row>` with prob columns (CSR-style) | Cache-friendly, O(rows) merge between H1 / H2 (parallel scan), easy lex-sort | Insertion in the middle costs O(rows); requires sort at end of construction | **chosen** |
| Indexed by global packed `int64` key | Compact, supports hashing | Same data, different key encoding; reserve for F4.x optimisation | secondary |

The chosen representation matches the R-reference (a `data.frame` with one
column per member sorted by `do.call(order, ...)`), so the row-by-row 1e-12
comparison in F2.6 is straightforward.

For lookups during composition (§8.5) the joint table is iterated linearly;
no hash map is required.

---

## 10. Bindings boundary

The binding wrapper (`src/rcpp_bindings.cpp` and, from F2.x onward,
`src/rcpp/*.cpp` if file growth justifies a split) is the *only* place
that:

- Translates R types to / from `core::*` types
  (`pedtools::ped` → `core::Pedigree`, named `numeric` → `core::Marker`,
  S3 `mutation` list → `core::MutationModel`, etc.).
- Calls `Rcpp::stop` on a `Result::error`.
- Optionally accepts an `arma::mat` for dense-matrix arguments via
  RcppArmadillo, then converts to `std::vector<std::vector<double>>` for
  core consumption (or, for performance-sensitive paths, the converse:
  core returns `arma::mat`-friendly contiguous storage that the wrapper
  zero-copies into Armadillo).

**F2.5 status:** the boundary now uses RcppArmadillo for the dense
return paths of `mutation_matrix_cpp` (K×K `arma::mat`) and
`cpt_marker_joint_cpp` (n_rows × n_members `arma::imat`), via two
private helpers `row_major_to_arma` / `states_flat_to_arma` defined at
the top of `src/rcpp_bindings.cpp`. The probability vectors (`P_H1`,
`P_H2`) stay as `Rcpp::NumericVector` because `arma::vec` wraps to a
1-column matrix in R (dim attribute), which is the wrong shape for the
R-side data.frame assignment in `R/cpt_marker_joint_cpp.R`. Bench
results: `mispitools_2_loop/benchmarks/F2.5_summary.md`.

The R-side public API (`R/marker_model.R`, `R/per_marker_kl.R` etc.) sits
*above* the binding layer and is the surface that user code calls. The
binding layer never gets called directly by user code — only by the R-side
functions, which guarantee precondition validation.

This three-layer split (user-R → boundary-R → bindings → core) is the
*Rcpp standard idiom*. The same pattern shows up later for WASM
(user-JS → boundary-JS → embind glue → core) and for the CLI
(`main.cpp` → core directly).

---

## 11. Cross-target portability

### 11.1 WASM (Emscripten, F10)

`core/` builds under `emcc -std=c++17 -O3` with no flags beyond
`-fno-exceptions` (WASM exceptions are expensive). The chosen error-handling
pattern (`Result<T>` with `std::optional` + `std::string`) is exception-free,
so `-fno-exceptions` does not break correctness.

Embind layer in `mispitools_2_loop/wasm/bindings.cpp` declares JS-visible
types `Pedigree`, `Marker`, `MutationModel`, `JointTable`, and the
top-level functions `cpt_marker_joint`, `per_marker_kl`, `per_marker_lr_dist`,
`evaluate_evidence`. The TS wrapper at `mispitools_2_loop/web/src/wasm/mispitools.ts`
gives users typed access.

Known pitfalls:

- `std::random_device` is not supported in WASM. The kernel doesn't sample,
  but tests must be aware. (Sampling is done from JS via `Math.random()` /
  `crypto.getRandomValues()` and passed in as a vector.)
- Long-running computations block the JS event loop. F10+ wraps
  `evaluate_evidence` in a Web Worker. The kernel sees nothing.
- `RcppArmadillo` is *not* available under emcc. Anything in `core/` that
  would have used Armadillo must use `std::vector` instead. This is why
  RcppArmadillo is *banned in core/* (§2).

### 11.2 Standalone gcc/clang (CLI, post-2.0)

The standalone smoke test in `mispitools_2_loop/scripts/test_core_standalone.sh`
already builds `core/` with `gcc -std=c++17 -Wall -Wextra -Wpedantic`. As
real functions land in F2.2+, the script will extend to call them and
assert the same 1e-12 oracle as the R tests.

A CLI binary (`mispitools-cli`) would link `core/` plus a slim
`cli/main.cpp` and accept JSON inputs. Out of scope for 2.0.

---

## 12. Test / oracle strategy

### 12.1 Inside-package tests

For each module M added by F2.x, three categories of tests live in
`tests/testthat/`:

- **C++ ↔ R cross-check** (`test-cpt-marker-joint-cpp.R`,
  `test-mutation-matrix-cpp.R`, ...). Calls the binding wrapper from R,
  compares element-by-element against the R reference engine to 1e-12.
  Mirrors the structure of F1's tests against `pedprobr` / `pedmut` / `fbnet`.
- **Numerical-correctness** against analytic cases (trio HWE, AAxAB,
  ABxAB, half-sib full enumeration). Tol 1e-12.
- **Boundary cases**: K=2, K=20, mutation `rate → 0` (continuity with
  `none`), missing-parent imputation, fully-typed pedigree, singleton-POI.

### 12.2 Standalone smoke tests

`mispitools_2_loop/scripts/test_core_standalone.sh` runs after each F2.x
landing. As functions move from placeholders to real, the script picks them
up. The contract: standalone build *and* test pass with the same numerical
output as the R build, to 1e-12.

### 12.3 External oracles (kept in Suggests)

| Oracle | Used in | Phase |
|---|---|---|
| `pedprobr::oneMarkerDistribution` | `cpt_marker_joint` w/ `none` | F2.6 (already in F1.3 for R-ref) |
| `pedmut::mutationMatrix` | mutation builders | F2.6 |
| `fbnet::getLocusCPT` | mutation builders + joint reshape | F2.6 (already in F1.5 for R-ref) |
| `Familias::FamiliasPosterior$LRperMarker` | full LR per marker | F4.5 |
| `Familias::FamiliasPosterior` with Custom matrix | asymmetric mutation | F5.1 |
| `forensIT::perMarkerKLs` | per-marker KL (after upstream bug fix) | F3.3 (kept as documented expected oracle) |
| `forrel::profileSim` Monte Carlo | LR distribution histogram | F4.5 |

`pedprobr`, `pedmut`, `fbnet`, `Familias`, `forensIT`, `forrel` all live in
Suggests with `skip_if_not_installed()`. After F9.2 the policy is kept
("oracle-only") indefinitely.

### 12.4 R reference engine — perennial oracle

`R/r_ref_cpt.R` + `R/r_ref_per_marker.R` are *not* removed when the C++
engine lands. They stay in the package as `@noRd` reference implementations
and continue to be the 1e-12 oracle in CI. They are the slowest possible
correct implementation; the speed differential to C++ defines our F7
performance gains.

---

## 13. Performance targets (re §"Performance targets" in ROADMAP)

| Case | R-ref wall | Target C++ wall | Speedup |
|---|---|---|---|
| Trio + 13 STRs + 3 NG, none | ~3 s | < 50 ms | ≥ 60× |
| First-cousin + 21 STRs + 5 NG, equal | ~30 s | < 150 ms | ≥ 200× |
| 3-gen + 23 STRs + 5 NG, stepwise | ~120 s | < 300 ms | ≥ 400× |
| Memory peak | ~500 MB (R-ref) | < 50 MB |  |

Hotspots to engineer for from F2.2 forward, in priority order:

1. Avoid materialising the row Cartesian during peeling — extend in
   place, dropping zero rows after each non-founder is added.
2. Allele lumping (F2.5) for STRs with ≥ 10 alleles. Typically reduces
   `n_alleles` to 4 or 5 with no loss for typed members.
3. Vectorised mutation-matrix application via Armadillo at the binding
   layer when `n_alleles ≥ 10`; small dense `std::vector<std::vector<double>>`
   in core for `n_alleles < 10`.
4. Per-marker independent loops trivially parallel (F7.1 OpenMP).
5. Convolution short-circuit when one operand is a delta (single-support
   distribution).

OpenMP (F7.1) lives in the binding layer or in a delimited `#ifdef
MISPI_HAS_OPENMP` block within `core/`. Either way the `core/` API stays
serial; parallelism is a wrapping concern.

---

## 14. Linkage (F5.2)

Two-marker joint when `theta < 0.5`:

```
P(g_A, g_B | H, theta) = sum over selectors S_A, S_B of:
    P(g_A | parental copies, S_A) ·
    P(g_B | parental copies, S_B) ·
    P(S_B | S_A, theta)
```

The selector CPT is a 2×2 matrix with entries `(1 - r, r ; r, 1 - r)` where
`r` is the recombination fraction (Haldane). Up to two linked markers are
supported in 2.0; arbitrary linkage groups are post-2.0.

API:

```cpp
Result<JointTable> linked_pair_joint(
    const Pedigree& p,
    const Marker& m_A,
    const Marker& m_B,
    const MutationModel& mu_A,
    const MutationModel& mu_B,
    double theta);
```

The joint table now has 2 × `n_relevant_members` columns (one genotype per
marker per member). KL and LR distribution apply unchanged to the
two-marker joint — they only see the row table.

---

## 15. Things explicitly out of scope for `core/` in 2.0

These are noted so the kernel does not gain accidental complexity:

- **Loopy pedigrees** (inbreeding, half-sib loops with shared ancestor).
  The R wrapper calls `pedtools::breakLoops` at the boundary; `core/` sees
  acyclic input only. Post-2.0: cutset conditioning lift from Familias
  (SCOUT_Familias_cpp.md "NO lift" section becomes "patrón a imitar").
- **Multi-pedigree posterior**. No prior over pedigrees, no posterior
  combination. The user supplies one pedigree per `mp_case`; LR computation
  is hypothesis-pair only.
- **Silent alleles**. F5+ feature.
- **Population substructure (Balding-Nichols θ)**. F6+. The lift (SCOUT
  Familias Lift 1) is ready to drop in when needed.
- **Mutation matrix stabilisation (PM, BA, MH transforms)**. F5.x. The
  pedmut lift (SCOUT_pedmut.md Lift 3) is ready.
- **Lumping with arbitrary partitions** (Kemeny-Snell `lumpedMatrix`). F2.5+.
- **Backwards compatibility with mispitools 1.x outputs**. None. 2.0 is a
  break, per ROADMAP §"Decisiones de scope" point 2.

---

## 16. Open questions deferred to later milestones

- Genotype encoding: stick with the upper-triangular `(a1*K - a1(a1-1)/2 + (a2-a1))`
  convention, or move to a packed `int64` (`a1 << 32 | a2`)? Decision in F2.2
  when the first real CPT lands. Current default: upper-triangular int32.
- `LrDist::has_pos_inf` / `has_neg_inf` flags vs sentinel rows: decide in F4.1.
- Whether the `MutationModel` POD should carry a precomputed K×K matrix
  (eager) or only the parameters (lazy, rebuilt on each call). Decision in
  F2.4 when caching becomes a measurable win.
- `evidence_combine_markov_se` exact formula — pending Egeland-Marsico 2026
  paper draft. Confirmed in F6.5.

---

## 17. References

- F1 R-reference engine: `R/r_ref_cpt.R`, `R/r_ref_per_marker.R`.
- F1.8 SCOUT notes:
  `mispitools_2_loop/scout/SCOUT_fbnet_internals.md`,
  `mispitools_2_loop/scout/SCOUT_Familias_cpp.md`,
  `mispitools_2_loop/scout/SCOUT_pedmut.md`.
- ROADMAP: `mispitools_2_loop/ROADMAP.md`.
- STATE table: `mispitools_2_loop/STATE.md`.
- Papers on the underlying math:
  - Brinkmann B. et al. (1998), *Mutation rate in human microsatellites*, AJHG 62(6).
  - Dawid A.P. (2002), *Properties of diagnostic data distributions*, Biometrics.
  - Kemeny J.G., Snell J.L. (1976), *Finite Markov Chains*, Springer, Ch. 6.
  - Mostad P., Egeland T. (2017), *Familias* manual.
  - Vigeland M.D. (2020), *Pedigree Analysis in R*, Academic Press.
  - Egeland T., Marsico F.L. (2026, in prep), *Belief dynamics in forensic genetics*.
