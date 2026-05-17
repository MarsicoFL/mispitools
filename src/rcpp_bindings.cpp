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
// engine cpt_marker_joint_R(). Wires kinds 0 (None), 1 (Equal),
// 2 (Stepwise), 4 (Asymmetric/Dawid, F5.1 — uses marker.freqs as afreq);
// kind 3 (Proportional) still errors out at build_mutation_matrix().
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

// ---------------------------------------------------------------------------
// F3.1 — cpp_per_marker_kl(): bidirectional KL + expected log10 LR over a
// pre-computed sparse joint (P_H1, P_H2).
//
// Takes the joint as two parallel numeric vectors (the columns produced by
// cpt_marker_joint_cpp / cpt_marker_joint_R). Returns a list with both KLs,
// both expectations, and the KLde-style absolute-continuity diagnostics.
// The R-side wrapper that consumes a marker_model directly arrives in F3.2.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_per_marker_kl(
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2) {
    if (p_h1.size() != p_h2.size()) {
        Rcpp::stop("cpp_per_marker_kl: p_h1 and p_h2 must have the same length.");
    }

    mc::JointTable jt;
    jt.n_members = 0;
    jt.n_genotypes = 0;
    jt.p_h1.assign(p_h1.begin(), p_h1.end());
    jt.p_h2.assign(p_h2.begin(), p_h2.end());

    auto r = mc::per_marker_kl(jt);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::PerMarkerKL& v = *r;

    return Rcpp::List::create(
        Rcpp::_["e_log10_lr_h1"]            = v.e_log10_lr_h1,
        Rcpp::_["e_log10_lr_h2"]            = v.e_log10_lr_h2,
        Rcpp::_["kl_h1_to_h2"]              = v.kl_h1_to_h2,
        Rcpp::_["kl_h2_to_h1"]              = v.kl_h2_to_h1,
        Rcpp::_["abs_cont_violations_h2"]   = static_cast<int>(v.abs_cont_violations_h2),
        Rcpp::_["abs_cont_violations_h1"]   = static_cast<int>(v.abs_cont_violations_h1),
        Rcpp::_["mass_violations_h2"]       = v.mass_violations_h2,
        Rcpp::_["mass_violations_h1"]       = v.mass_violations_h1
    );
}

// ---------------------------------------------------------------------------
// F4.1 — cpp_per_marker_lr_dist(): sparse per-marker LR distribution over a
// pre-computed joint (P_H1, P_H2).
//
// Same input contract as cpp_per_marker_kl (the columns produced by
// cpt_marker_joint_cpp / cpt_marker_joint_R). Returns the sparse-sorted
// (log10_lr, p_h1, p_h2) atoms plus the ±Inf flags. `aggregate = TRUE`
// (default) collapses equal-log10_lr atoms, bit-for-bit with the R
// reference per_marker_lr_dist_R(). The model-aware R wrapper arrives in
// F4.4 (lr_distribution()).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_per_marker_lr_dist(
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2,
        bool aggregate = true) {
    if (p_h1.size() != p_h2.size()) {
        Rcpp::stop("cpp_per_marker_lr_dist: p_h1 and p_h2 must have the same length.");
    }

    mc::JointTable jt;
    jt.n_members = 0;
    jt.n_genotypes = 0;
    jt.p_h1.assign(p_h1.begin(), p_h1.end());
    jt.p_h2.assign(p_h2.begin(), p_h2.end());

    auto r = mc::per_marker_lr_dist(jt, aggregate);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::LrDist& d = *r;

    return Rcpp::List::create(
        Rcpp::_["log10_lr"]    = Rcpp::NumericVector(d.log10_lr.begin(), d.log10_lr.end()),
        Rcpp::_["p_h1"]        = Rcpp::NumericVector(d.p_h1.begin(), d.p_h1.end()),
        Rcpp::_["p_h2"]        = Rcpp::NumericVector(d.p_h2.begin(), d.p_h2.end()),
        Rcpp::_["has_pos_inf"] = d.has_pos_inf,
        Rcpp::_["has_neg_inf"] = d.has_neg_inf
    );
}

// ---------------------------------------------------------------------------
// F4.2 — cpp_lr_dist_compose(): convolution of independent per-feature LR
// distributions.
//
// `dists` is a list of lists, each with numeric `log10_lr`, `p_h1`, `p_h2`
// columns (the per-feature output of cpp_per_marker_lr_dist). `method` is
// "exact" (default) or "grid". Returns the composed sparse-sorted
// distribution + the ±Inf flags. The model-aware R wrapper that builds the
// per-feature list from a model list arrives in F4.4 (lr_distribution()).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_lr_dist_compose(
        Rcpp::List dists,
        std::string method = "exact",
        double merge_tol = 0.0,
        int grid_points = 512) {
    mc::LrDistComposeOptions opts;
    if (method == "exact") {
        opts.method = mc::ComposeMethod::Exact;
    } else if (method == "grid") {
        opts.method = mc::ComposeMethod::Grid;
    } else {
        Rcpp::stop("cpp_lr_dist_compose: method must be 'exact' or 'grid'.");
    }
    opts.merge_tol = merge_tol;
    opts.grid_points = grid_points;

    std::vector<mc::LrDist> per_feature;
    per_feature.reserve(static_cast<std::size_t>(dists.size()));
    for (R_xlen_t i = 0; i < dists.size(); ++i) {
        Rcpp::List d = dists[i];
        Rcpp::NumericVector lr  = d["log10_lr"];
        Rcpp::NumericVector p1  = d["p_h1"];
        Rcpp::NumericVector p2  = d["p_h2"];
        mc::LrDist x;
        x.log10_lr.assign(lr.begin(), lr.end());
        x.p_h1.assign(p1.begin(), p1.end());
        x.p_h2.assign(p2.begin(), p2.end());
        per_feature.push_back(std::move(x));
    }

    auto r = mc::lr_dist_compose(per_feature, opts);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::LrDist& d = *r;

    return Rcpp::List::create(
        Rcpp::_["log10_lr"]    = Rcpp::NumericVector(d.log10_lr.begin(), d.log10_lr.end()),
        Rcpp::_["p_h1"]        = Rcpp::NumericVector(d.p_h1.begin(), d.p_h1.end()),
        Rcpp::_["p_h2"]        = Rcpp::NumericVector(d.p_h2.begin(), d.p_h2.end()),
        Rcpp::_["has_pos_inf"] = d.has_pos_inf,
        Rcpp::_["has_neg_inf"] = d.has_neg_inf
    );
}

// ---------------------------------------------------------------------------
// F4.3 — decision-theoretic primitives over a sparse LR distribution.
//
// Inputs are the (log10_lr, p_h1, p_h2) columns of cpp_per_marker_lr_dist
// / cpp_lr_dist_compose. The model-aware S3 layer (as_lr_dist + summary /
// plot methods) lives in R/lr_dist_s3.R; F4.4 adds lr_distribution().
// ---------------------------------------------------------------------------

namespace {

inline mc::LrDist lr_dist_from_cols(const Rcpp::NumericVector& lr,
                                    const Rcpp::NumericVector& p1,
                                    const Rcpp::NumericVector& p2) {
    mc::LrDist d;
    d.log10_lr.assign(lr.begin(), lr.end());
    d.p_h1.assign(p1.begin(), p1.end());
    d.p_h2.assign(p2.begin(), p2.end());
    return d;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_lr_dist_summary(
        Rcpp::NumericVector log10_lr,
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2) {
    auto r = mc::lr_dist_summary(lr_dist_from_cols(log10_lr, p_h1, p_h2));
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::LrDistSummary& s = *r;
    return Rcpp::List::create(
        Rcpp::_["mean_h1"]     = s.mean_h1,
        Rcpp::_["mean_h2"]     = s.mean_h2,
        Rcpp::_["var_h1"]      = s.var_h1,
        Rcpp::_["var_h2"]      = s.var_h2,
        Rcpp::_["sd_h1"]       = s.sd_h1,
        Rcpp::_["sd_h2"]       = s.sd_h2,
        Rcpp::_["mass_h1"]     = s.mass_h1,
        Rcpp::_["mass_h2"]     = s.mass_h2,
        Rcpp::_["has_pos_inf"] = s.has_pos_inf,
        Rcpp::_["has_neg_inf"] = s.has_neg_inf
    );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_lr_dist_quantile(
        Rcpp::NumericVector log10_lr,
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2,
        Rcpp::NumericVector probs,
        bool under_h1 = true) {
    std::vector<double> pv(probs.begin(), probs.end());
    auto r = mc::lr_dist_quantile(
        lr_dist_from_cols(log10_lr, p_h1, p_h2), pv, under_h1);
    if (!r.ok()) Rcpp::stop(r.error);
    const std::vector<double>& q = *r;
    return Rcpp::NumericVector(q.begin(), q.end());
}

// [[Rcpp::export]]
Rcpp::List cpp_lr_dist_decision_rates(
        Rcpp::NumericVector log10_lr,
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2,
        double threshold) {
    auto r = mc::decision_rates(
        lr_dist_from_cols(log10_lr, p_h1, p_h2), threshold);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::DecisionRates& v = *r;
    return Rcpp::List::create(
        Rcpp::_["threshold"] = v.threshold,
        Rcpp::_["fpr"]       = v.fpr,
        Rcpp::_["fnr"]       = v.fnr,
        Rcpp::_["tpr"]       = v.tpr,
        Rcpp::_["tnr"]       = v.tnr,
        Rcpp::_["mcc"]       = v.mcc
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_lr_dist_roc(
        Rcpp::NumericVector log10_lr,
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2) {
    auto r = mc::roc_curve(lr_dist_from_cols(log10_lr, p_h1, p_h2));
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::RocCurve& c = *r;
    return Rcpp::List::create(
        Rcpp::_["threshold"] = Rcpp::NumericVector(c.threshold.begin(), c.threshold.end()),
        Rcpp::_["fpr"]       = Rcpp::NumericVector(c.fpr.begin(), c.fpr.end()),
        Rcpp::_["tpr"]       = Rcpp::NumericVector(c.tpr.begin(), c.tpr.end()),
        Rcpp::_["fnr"]       = Rcpp::NumericVector(c.fnr.begin(), c.fnr.end()),
        Rcpp::_["tnr"]       = Rcpp::NumericVector(c.tnr.begin(), c.tnr.end()),
        Rcpp::_["auc"]       = c.auc
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_lr_dist_choose_threshold(
        Rcpp::NumericVector log10_lr,
        Rcpp::NumericVector p_h1,
        Rcpp::NumericVector p_h2,
        double weight = 10.0) {
    auto r = mc::choose_threshold_weighted(
        lr_dist_from_cols(log10_lr, p_h1, p_h2), weight);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::ThresholdChoice& v = *r;
    return Rcpp::List::create(
        Rcpp::_["threshold"] = v.threshold,
        Rcpp::_["fpr"]       = v.fpr,
        Rcpp::_["fnr"]       = v.fnr,
        Rcpp::_["distance"]  = v.distance
    );
}

// ---------------------------------------------------------------------------
// F3.4 — cpp_per_marker_kl_batch(): N-marker batch over one shared pedigree.
//
// All markers share the topology (`father`, `mother`, `poi`). For each
// marker the binding flattens freqs + numeric_labels + mutation parameters
// into POD vectors and the core loops once, caching mutation matrices
// across markers with identical (kind, K, rate, range, labels) keys.
// Returns the per-marker KL columns + cache diagnostics.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_per_marker_kl_batch(
        Rcpp::IntegerVector father,
        Rcpp::IntegerVector mother,
        int poi,
        Rcpp::List freqs_list,
        Rcpp::IntegerVector mutation_kind,
        Rcpp::NumericVector mutation_rate,
        Rcpp::NumericVector mutation_range,
        Rcpp::List numeric_labels_list) {
    if (father.size() != mother.size()) {
        Rcpp::stop("cpp_per_marker_kl_batch: father and mother must be the same length.");
    }
    const R_xlen_t N = freqs_list.size();
    if (mutation_kind.size() != N || mutation_rate.size() != N
            || mutation_range.size() != N
            || numeric_labels_list.size() != N) {
        Rcpp::stop("cpp_per_marker_kl_batch: per-marker vectors must all have the same length.");
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

    std::vector<mc::Marker> markers;
    std::vector<mc::MutationModel> mutations;
    markers.reserve(static_cast<std::size_t>(N));
    mutations.reserve(static_cast<std::size_t>(N));

    for (R_xlen_t i = 0; i < N; ++i) {
        Rcpp::NumericVector f = freqs_list[i];
        Rcpp::NumericVector lab = numeric_labels_list[i];
        const int kind_i = mutation_kind[i];
        if (kind_i < 0 || kind_i > 4) {
            Rcpp::stop("cpp_per_marker_kl_batch: mutation_kind out of range [0, 4] at marker %d.",
                       static_cast<int>(i + 1));
        }

        mc::Marker m;
        m.id = "";
        m.n_alleles = static_cast<mc::AlleleIndex>(f.size());
        m.freqs.assign(f.begin(), f.end());
        m.numeric_labels.assign(lab.begin(), lab.end());
        markers.push_back(std::move(m));

        mc::MutationModel mut;
        mut.kind  = static_cast<mc::MutationKind>(kind_i);
        mut.rate  = mutation_rate[i];
        mut.range = mutation_range[i];
        mut.rate2 = 0.0;
        mut.bias  = 0.5;
        mutations.push_back(mut);
    }

    auto r = mc::per_marker_kl_batch(ped, markers, mutations);
    if (!r.ok()) Rcpp::stop(r.error);
    const mc::PerMarkerKLBatch& batch = *r;

    Rcpp::NumericVector e_h1(N), e_h2(N), kl_12(N), kl_21(N);
    Rcpp::IntegerVector vio_h1(N), vio_h2(N);
    Rcpp::NumericVector mass_h1(N), mass_h2(N);

    for (R_xlen_t i = 0; i < N; ++i) {
        const mc::PerMarkerKL& v = batch.entries[static_cast<std::size_t>(i)];
        e_h1[i]   = v.e_log10_lr_h1;
        e_h2[i]   = v.e_log10_lr_h2;
        kl_12[i]  = v.kl_h1_to_h2;
        kl_21[i]  = v.kl_h2_to_h1;
        vio_h1[i] = static_cast<int>(v.abs_cont_violations_h1);
        vio_h2[i] = static_cast<int>(v.abs_cont_violations_h2);
        mass_h1[i] = v.mass_violations_h1;
        mass_h2[i] = v.mass_violations_h2;
    }

    return Rcpp::List::create(
        Rcpp::_["e_log10_lr_h1"]          = e_h1,
        Rcpp::_["e_log10_lr_h2"]          = e_h2,
        Rcpp::_["kl_h1_to_h2"]            = kl_12,
        Rcpp::_["kl_h2_to_h1"]            = kl_21,
        Rcpp::_["abs_cont_violations_h1"] = vio_h1,
        Rcpp::_["abs_cont_violations_h2"] = vio_h2,
        Rcpp::_["mass_violations_h1"]     = mass_h1,
        Rcpp::_["mass_violations_h2"]     = mass_h2,
        Rcpp::_["cache_hits"]   = static_cast<int>(batch.mutation_matrix_cache_hits),
        Rcpp::_["cache_misses"] = static_cast<int>(batch.mutation_matrix_cache_misses)
    );
}

// [[Rcpp::export]]
arma::mat mutation_matrix_cpp(
        int K,
        int mutation_kind = 0,
        double mutation_rate = 0.0,
        double mutation_range = 0.0,
        Rcpp::NumericVector numeric_labels = Rcpp::NumericVector::create(),
        Rcpp::NumericVector afreq = Rcpp::NumericVector::create()) {
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
    std::vector<double> p(afreq.begin(), afreq.end());
    auto r = mc::build_mutation_matrix(
        mut, static_cast<mc::AlleleIndex>(K), labels, p);
    if (!r.ok()) Rcpp::stop(r.error);

    // F2.5 — row-major (core) → column-major (R) via Armadillo's
    // vectorised transpose. Replaces the previous O(K^2) element-wise
    // assignment loop into Rcpp::NumericMatrix.
    return row_major_to_arma(*r, static_cast<arma::uword>(K));
}

// ---------------------------------------------------------------------------
// F5.1 — dawid_max_rate_cpp(): largest well-defined `rate` for the
// asymmetric (Dawid 2002) model given allele frequencies and range.
// Mirrors pedmut::maxRate() (UW bound). Used by the verifier to confirm
// the undefined-model cap matches the oracle within tolerance.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double dawid_max_rate_cpp(Rcpp::NumericVector afreq, double range) {
    std::vector<double> p(afreq.begin(), afreq.end());
    auto r = mc::dawid_max_rate(p, range);
    if (!r.ok()) Rcpp::stop(r.error);
    return *r;
}
