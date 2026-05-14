// F0.5 — Rcpp bindings over the pure-C++ core in src/core/.
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

#include <Rcpp.h>

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
// F2.2 — cpt_marker_joint_cpp(): joint genotype CPT under H1 / H2.
//
// The binding flattens the marker_model + pedtools::ped object into POD
// vectors on the R side (see R/cpt_marker_joint_cpp.R) and the core
// computes the sparse joint table. The R-side wrapper rebuilds the
// data.frame view with labelled genotype strings to match the R-reference
// engine cpt_marker_joint_R(). Mutation kinds other than 0 (None) return
// an error in F2.2; Equal/Stepwise arrive in F2.4.
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

    const int n_rows = static_cast<int>(jt.p_h1.size());
    const int n_mem = ped.n_members;

    // 1-based genotype indices on the R side (matches the R-reference engine).
    Rcpp::IntegerMatrix states(n_rows, n_mem);
    for (int r_ = 0; r_ < n_rows; ++r_) {
        for (int i = 0; i < n_mem; ++i) {
            states(r_, i) = jt.states_flat[
                static_cast<std::size_t>(r_) * n_mem + i] + 1;
        }
    }
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
// F2.3 — mutation_matrix_cpp(): K x K mutation matrix for None/Equal/Stepwise.
//
// Dispatches to the pure-core builders in mutation_models.cpp. Returns a
// K x K NumericMatrix in row-major-equivalent layout (R matrix populated
// from M[i, j] = mat[i * K + j]). Mutation kinds outside {0, 1, 2} raise
// an R error to surface the boundary clearly; Asymmetric (4) arrives in
// F5.1.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix mutation_matrix_cpp(
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
    const std::vector<double>& flat = *r;

    Rcpp::NumericMatrix out(K, K);
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            out(i, j) = flat[static_cast<std::size_t>(i) * K + j];
        }
    }
    return out;
}
