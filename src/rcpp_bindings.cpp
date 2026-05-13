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
