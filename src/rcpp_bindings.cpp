// F0.5 / F2.5 — Rcpp(+Armadillo) bindings over the pure-C++ core in src/core/.
//
// Layered architecture for mispitools 2.0:
//
//   src/core/             pure C++17, no Rcpp/R headers, portable to WASM.
//   src/rcpp_bindings.cpp this file — thin wrappers crossing the R boundary.
//
// Note on location: the ROADMAP sketches the binding layer under
// src/rcpp/ but Rcpp::compileAttributes() only scans the top level of
// src/ (no recursion), and devtools::load_all() invokes it implicitly
// at every build. Keeping the wrappers at top-level src/ lets the
// generated src/RcppExports.cpp and R/RcppExports.R stay in sync
// automatically. The pure C++ engine remains in src/core/ as planned.
//
// F2.5 — RcppArmadillo is enabled here (not in src/core/, which stays
// portable to Emscripten — see DESIGN.md §2 / §11.1). The boundary uses
// `arma::mat` / `arma::imat` for the dense outputs of mutation_matrix_cpp
// and cpt_marker_joint_cpp so the row-major → column-major transpose runs
// through Armadillo's vectorised copy + `Rcpp::wrap()` zero-copy
// conversion, instead of an element-wise R-side `out(i, j) = ...` loop.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "core/pedigree.h"
#include "core/marker.h"
#include "core/mutation_models.h"
#include "core/cpt_engine.h"
#include "core/kl_engine.h"
#include "core/lr_dist.h"
#include "core/nongenetic_lr.h"
#include "core/evidence_combine.h"
#include "core/concentration.h"
#include "core/decision.h"
#include "core/linkage.h"

namespace mc = mispitools::core;

namespace {

// Row-major std::vector<double> (length K*K, M[i,j] = flat[i*K + j])
// → arma::mat (K x K). Armadillo is column-major, so we read the buffer
// into a column-major view (advisory_ok = false, copy_aux_mem = true to
// own the memory) and then transpose. The transpose is a contiguous
// vectorised copy under Armadillo's standard kernel, faster than the
// cell-wise loop the F2.3 binding used.
inline arma::mat row_major_to_arma(
        const std::vector<double>& flat,
        arma::uword K) {
    if (K == 0) return arma::mat();
    arma::mat view(const_cast<double*>(flat.data()), K, K,
                   /*copy_aux_mem=*/false, /*strict=*/true);
    return view.t();
}

// Row-major std::vector<GenotypeIndex> states_flat (n_rows * n_members)
// → arma::imat (n_rows x n_members), 1-based to match R-side genotype
// indices in `cpt_marker_joint_R`. Same trick: column-major view +
// transpose, then unary `+ 1`.
inline arma::imat states_flat_to_arma(
        const std::vector<mc::GenotypeIndex>& flat,
        arma::uword n_rows,
        arma::uword n_members) {
    if (n_rows == 0 || n_members == 0) return arma::imat();
    // GenotypeIndex is std::int32_t; arma::imat element is sword
    // (typically long long). Copy through a transposed view rather
    // than reinterpret, since the element widths differ.
    arma::imat out(n_rows, n_members);
    // Fill column-by-column; each column reads strided from the
    // row-major source (member i runs at offset i, stride n_members).
    for (arma::uword i = 0; i < n_members; ++i) {
        arma::sword* col = out.colptr(i);
        const mc::GenotypeIndex* src = flat.data() + i;
        for (arma::uword r = 0; r < n_rows; ++r) {
            col[r] = static_cast<arma::sword>(src[r * n_members]) + 1;
        }
    }
    return out;
}

}  // namespace

// [[Rcpp::export]]
int cpp_pedigree_placeholder(int x)         { return mc::pedigree_placeholder(x); }

// [[Rcpp::export]]
int cpp_marker_placeholder(int x)           { return mc::marker_placeholder(x); }

// [[Rcpp::export]]
int cpp_mutation_models_placeholder(int x)  { return mc::mutation_models_placeholder(x); }

// [[Rcpp::export]]
int cpp_cpt_engine_placeholder(int x)       { return mc::cpt_engine_placeholder(x); }

// [[Rcpp::export]]
int cpp_kl_engine_placeholder(int x)        { return mc::kl_engine_placeholder(x); }

// [[Rcpp::export]]
int cpp_lr_dist_placeholder(int x)          { return mc::lr_dist_placeholder(x); }

// [[Rcpp::export]]
int cpp_nongenetic_lr_placeholder(int x)    { return mc::nongenetic_lr_placeholder(x); }

// [[Rcpp::export]]
int cpp_evidence_combine_placeholder(int x) { return mc::evidence_combine_placeholder(x); }

// [[Rcpp::export]]
int cpp_concentration_placeholder(int x)    { return mc::concentration_placeholder(x); }

// [[Rcpp::export]]
int cpp_decision_placeholder(int x)         { return mc::decision_placeholder(x); }

// [[Rcpp::export]]
int cpp_linkage_placeholder(int x)          { return mc::linkage_placeholder(x); }

// ---------------------------------------------------------------------------
// F2.2 / F2.4 — cpt_marker_joint_cpp(): joint genotype CPT under H1 / H2.
//
// The binding flattens the marker_model + pedtools::ped object into POD
// vectors on the R side (see R/cpt_marker_joint_cpp.R) and the core
// computes the sparse joint table. The R-side wrapper rebuilds the
// data.frame view with labelled genotype strings to match the R-reference
// engine cpt_marker_joint_R(). F2.4 wires kinds 0 (None), 1 (Equal),
// 2 (Stepwise); kinds 3 (Proportional) and 4 (Asymmetric) error out at
// build_mutation_matrix() until F5.1.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpt_marker_joint_cpp(
        Rcpp::IntegerVector father,
        Rcpp::IntegerVector mother,
        int poi,
        Rcpp::NumericVector freqs,
        int mutation_kind = 0,
        double mutation_rate = 0.0,
        double mutation_range = 0.0,
        double mutation_rate2 = 0.0,
        double mutation_bias = 0.5,
        Rcpp::NumericVector numeric_labels = Rcpp::NumericVector::create()) {
    if (father.size() != mother.size()) {
        Rcpp::stop("cpt_marker_joint_cpp: father and mother must be the same length.");
    }
    if (mutation_kind < 0 || mutation_kind > 4) {
        Rcpp::stop("cpt_marker_joint_cpp: mutation_kind out of range [0, 4].");
    }

    mc::Pedigree ped;
    ped.n_members = static_cast<mc::MemberIndex>(father.size());
    ped.father.assign(static_cast<std::size_t>(ped.n_members), mc::kNoParent);
    ped.mother.assign(static_cast<std::size_t>(ped.n_members), mc::kNoParent);
    for (int i = 0; i < ped.n_members; ++i) {
        ped.father[static_cast<std::size_t>(i)] = father[i];
        ped.mother[static_cast<std::size_t>(i)] = mother[i];
    }
    ped.poi = poi;

    mc::Marker marker;
    marker.id = "";
    marker.n_alleles = static_cast<mc::AlleleIndex>(freqs.size());
    marker.freqs.assign(freqs.begin(), freqs.end());
    marker.numeric_labels.assign(numeric_labels.begin(), numeric_labels.end());

    mc::MutationModel mut;
    mut.kind = static_cast<mc::MutationKind>(mutation_kind);
    mut.rate = mutation_rate;
    mut.range = mutation_range;
    mut.rate2 = mutation_rate2;
    mut.bias = mutation_bias;

    auto r = mc::cpt_marker_joint(ped, marker, mut);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::JointTable& jt = *r;

    const arma::uword n_rows = jt.p_h1.size();
    const arma::uword n_mem = static_cast<arma::uword>(ped.n_members);

    // F2.5 — row-major (core) → column-major (R) handled by an Armadillo
    // strided copy inside states_flat_to_arma. arma::imat → R `integer`
    // matrix via Rcpp::wrap (zero-copy of Armadillo's storage). The
    // probability vectors stay as plain Rcpp::NumericVector — arma::vec
    // wraps to a 1-column matrix (dim attribute) which is the wrong
    // shape for the R-side `out$P_H1 <- res$P_H1` assignment in
    // R/cpt_marker_joint_cpp.R.
    arma::imat states = states_flat_to_arma(jt.states_flat, n_rows, n_mem);
    Rcpp::NumericVector p_h1(jt.p_h1.begin(), jt.p_h1.end());
    Rcpp::NumericVector p_h2(jt.p_h2.begin(), jt.p_h2.end());

    return Rcpp::List::create(
        Rcpp::_["states"] = states,
        Rcpp::_["P_H1"]   = p_h1,
        Rcpp::_["P_H2"]   = p_h2,
        Rcpp::_["n_genotypes"] = static_cast<int>(jt.n_genotypes)
    );
}

// ---------------------------------------------------------------------------
// F2.3 / F2.5 — mutation_matrix_cpp(): K x K mutation matrix for
// None/Equal/Stepwise.
//
// Dispatches to the pure-core builders in mutation_models.cpp. F2.5
// returns an arma::mat carrying the same `M[i, j]` semantics as the
// previous Rcpp::NumericMatrix; the row-major → column-major conversion
// happens once inside `row_major_to_arma()` via an Armadillo transpose.
// Mutation kinds outside {0, 1, 2} raise an R error to surface the
// boundary clearly; Asymmetric (4) arrives in F5.1.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat mutation_matrix_cpp(
        int K,
        int mutation_kind = 0,
        double mutation_rate = 0.0,
        double mutation_range = 0.0,
        Rcpp::NumericVector numeric_labels = Rcpp::NumericVector::create()) {
    if (K <= 0) {
        Rcpp::stop("mutation_matrix_cpp: K must be positive.");
    }
    if (mutation_kind < 0 || mutation_kind > 4) {
        Rcpp::stop("mutation_matrix_cpp: mutation_kind out of range [0, 4].");
    }

    mc::MutationModel mut;
    mut.kind = static_cast<mc::MutationKind>(mutation_kind);
    mut.rate = mutation_rate;
    mut.range = mutation_range;

    std::vector<double> labels(numeric_labels.begin(), numeric_labels.end());
    auto r = mc::build_mutation_matrix(
        mut, static_cast<mc::AlleleIndex>(K), labels);
    if (!r.ok()) Rcpp::stop(r.error);

    // F2.5 — row-major (core) → column-major (R) via Armadillo's
    // vectorised transpose. Replaces the previous O(K^2) element-wise
    // assignment loop into Rcpp::NumericMatrix.
    return row_major_to_arma(*r, static_cast<arma::uword>(K));
}
