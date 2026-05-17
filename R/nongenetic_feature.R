#' Non-genetic feature model (S3)
#'
#' @description
#' Constructor for a single non-genetic forensic feature (sex, age, region,
#' hair colour, eye colour, pigmentation, birth date, or a user-defined
#' `"custom"` feature). It is the non-genetic counterpart of
#' [marker_model()]: it bundles the recorded observation, the
#' reference/observation sub-model, the population reference, and the
#' observation-error model into a single object that the per-feature engines
#' (arriving in milestones F6.3--F6.5) consume.
#'
#' Like [marker_model()], this constructor only validates and stores its
#' inputs; no conditional probability table is built here. The point of F6.1
#' is to fix the structural contract so that the genetic and non-genetic
#' paths are symmetric: every non-genetic feature reduces to a CPT under
#' `H1` (missing person, with observation error) and a CPT under `H2`
#' (population marginal), then `LR = P(D | H1) / P(D | H2)` --- exactly the
#' shape of the per-marker genetic path.
#'
#' @section Feature classes:
#' Each `type` maps to one of three internal *feature classes* that
#' determine how `db_or_freqs` and `error` are interpreted:
#' \describe{
#'   \item{`categorical`}{`sex`, `region`, `hair`, `eyes`, `pigmentation`.
#'     `db_or_freqs` is a named numeric vector of population category
#'     frequencies (the `H2` marginal); `error` is either a scalar symmetric
#'     misclassification rate or a full row-stochastic confusion matrix
#'     `E` with `E[true, observed]`. Discrete support, so the categorical
#'     KL is identical to the genetic engine.}
#'   \item{`continuous`}{`age`. The reference is either `"uniform"` over a
#'     numeric `range` or `"empirical"` from a sample passed in
#'     `db_or_freqs`; `error` is a scalar mis-binning rate.}
#'   \item{`date`}{`birthdate`. A Dirichlet model over signed
#'     declared-minus-actual day discrepancy bins (`cuts`); `error` is the
#'     Dirichlet `alpha` vector. `search = "open"` uses a uniform `H2`;
#'     `search = "closed"` uses database bin frequencies.}
#' }
#' `type = "custom"` requires `model` to declare its `class` (one of the
#' three above); validation then follows that class.
#'
#' @param type Character scalar. One of `"sex"`, `"age"`, `"region"`,
#'   `"hair"`, `"eyes"`, `"pigmentation"`, `"birthdate"`, `"custom"`.
#' @param observed The recorded observation for the unidentified person.
#'   For categorical features a length-one value matching a category of
#'   `db_or_freqs` (numeric labels are matched as characters). For `age`
#'   a single finite numeric. For `birthdate` either a single `Date` / a
#'   `"YYYY-MM-DD"` string, or a single finite numeric day discrepancy.
#' @param model Either `NULL` (a type-appropriate default is filled in) or
#'   a named list describing the reference/observation sub-model. For
#'   categorical features the field `reference` is one of `"marginal"`
#'   (use `db_or_freqs` as `H2`, the default) or `"uniform"`. For `age`,
#'   `reference` is `"uniform"` (with a length-2 increasing numeric
#'   `range`, default `c(1, 80)`) or `"empirical"`. For `birthdate`,
#'   `search` is `"open"` (default) or `"closed"`, with a strictly
#'   increasing numeric `cuts` vector (default
#'   `c(-120, -30, 30, 120, 240, 360)`). For `"custom"`, `model` must be a
#'   named list with a character `class` field in
#'   `c("categorical", "continuous", "date")`.
#' @param db_or_freqs The population reference. Categorical: a named numeric
#'   vector (>= 2 categories, entries in `[0, 1]`, summing to 1 within
#'   `tol`). Continuous: `NULL` when `reference = "uniform"`, otherwise a
#'   numeric sample (length >= 2, finite). Date: `NULL` for an open search,
#'   otherwise a non-negative numeric vector of bin frequencies of length
#'   `length(cuts) + 1` or a `data.frame` of declared dates.
#' @param error The observation-error model under `H1`. Categorical: a
#'   scalar in `[0, 1)` (symmetric misclassification) or a square
#'   row-stochastic numeric matrix whose dimension equals the number of
#'   categories. The matrix is positional: row `i` / column `j` are the
#'   `i`-th / `j`-th category in the order of `db_or_freqs` (any dimnames
#'   are informational only). Continuous: a scalar in `[0, 1)`. Date: a
#'   numeric
#'   Dirichlet `alpha` vector, all strictly positive, of length
#'   `length(cuts) + 1`. Defaults are class-appropriate except for
#'   `"custom"`, where `error` is required.
#' @param tol Numeric tolerance for the sum-to-one and row-stochastic
#'   checks. Defaults to `1e-6`.
#'
#' @return An object of class `"nongenetic_feature"`: a list with
#'   components `type`, `feature_class`, `observed`, `model`, `error`,
#'   `db_or_freqs`, and `categories` (the category labels for categorical
#'   features, `NULL` otherwise).
#'
#' @seealso [marker_model()] for the genetic counterpart; the legacy
#'   non-genetic functions [lr_sex()], [lr_age()], [lr_hair_color()],
#'   [lr_birthdate()], [lr_pigmentation()] whose behaviour this framework
#'   generalises and which serve as its regression oracles.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Categorical: biological sex, missing person female, 5% error
#' f_sex <- nongenetic_feature(
#'   type = "sex",
#'   observed = "F",
#'   db_or_freqs = c(F = 0.5, M = 0.5),
#'   error = 0.05
#' )
#' print(f_sex)
#'
#' # Categorical hair colour with a full confusion matrix
#' E <- error_matrix_hair()
#' f_hair <- nongenetic_feature(
#'   type = "hair",
#'   observed = 1,
#'   db_or_freqs = c("1" = 0.3, "2" = 0.2, "3" = 0.25, "4" = 0.15, "5" = 0.1),
#'   error = E
#' )
#'
#' # Continuous: age, uniform reference over [1, 80]
#' f_age <- nongenetic_feature(
#'   type = "age",
#'   observed = 42,
#'   model = list(reference = "uniform", range = c(1, 80)),
#'   error = 0.05
#' )
#'
#' # Date: birth-date discrepancy, open search, default Dirichlet
#' f_bd <- nongenetic_feature(
#'   type = "birthdate",
#'   observed = 45,
#'   error = c(1, 4, 60, 11, 6, 4, 4)
#' )
nongenetic_feature <- function(type,
                               observed,
                               model = NULL,
                               db_or_freqs = NULL,
                               error = NULL,
                               tol = 1e-6) {
  type <- validate_ng_type(type)
  fclass <- ng_feature_class(type, model)
  model <- validate_ng_model(model, type, fclass)
  db_or_freqs <- validate_ng_db(db_or_freqs, fclass, model, tol)
  cats <- if (fclass == "categorical") names(db_or_freqs) else NULL
  error <- validate_ng_error(error, fclass, cats, model, tol)
  observed <- validate_ng_observed(observed, fclass, cats)

  structure(
    list(
      type = type,
      feature_class = fclass,
      observed = observed,
      model = model,
      error = error,
      db_or_freqs = db_or_freqs,
      categories = cats
    ),
    class = "nongenetic_feature"
  )
}

#' @noRd
ng_types <- function() {
  c("sex", "age", "region", "hair", "eyes", "pigmentation",
    "birthdate", "custom")
}

#' @noRd
ng_classes <- function() {
  c("categorical", "continuous", "date")
}

#' @noRd
validate_ng_type <- function(type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !nzchar(type)) {
    stop("`type` must be a single non-empty character string.",
         call. = FALSE)
  }
  allowed <- ng_types()
  if (!type %in% allowed) {
    stop("`type` must be one of: ",
         paste(shQuote(allowed), collapse = ", "), ".", call. = FALSE)
  }
  type
}

#' @noRd
ng_feature_class <- function(type, model) {
  if (type == "custom") {
    if (is.null(model)) {
      stop("For `type = \"custom\"`, `model` must be supplied (a named ",
           "list with a `class` field).", call. = FALSE)
    }
    if (!is.list(model) || is.null(names(model)) ||
        is.null(model$class)) {
      stop("For `type = \"custom\"`, `model` must be a named list with a ",
           "`class` field.", call. = FALSE)
    }
    cls <- model$class
    if (!is.character(cls) || length(cls) != 1L || is.na(cls) ||
        !cls %in% ng_classes()) {
      stop("`model$class` must be one of: ",
           paste(shQuote(ng_classes()), collapse = ", "), ".",
           call. = FALSE)
    }
    return(cls)
  }
  switch(type,
    sex = "categorical",
    region = "categorical",
    hair = "categorical",
    eyes = "categorical",
    pigmentation = "categorical",
    age = "continuous",
    birthdate = "date"
  )
}

#' @noRd
default_ng_model <- function(type, fclass) {
  if (type == "custom") {
    # Class already validated; supply a class-appropriate skeleton.
    if (fclass == "categorical") {
      return(list(class = "categorical", reference = "marginal"))
    }
    if (fclass == "continuous") {
      return(list(class = "continuous", reference = "uniform"))
    }
    return(list(class = "date", search = "open",
                cuts = c(-120, -30, 30, 120, 240, 360)))
  }
  switch(type,
    sex = ,
    region = ,
    hair = ,
    eyes = ,
    pigmentation = list(reference = "marginal"),
    age = list(reference = "uniform", range = c(1, 80)),
    birthdate = list(search = "open",
                     cuts = c(-120, -30, 30, 120, 240, 360))
  )
}

#' @noRd
validate_ng_model <- function(model, type, fclass) {
  if (is.null(model)) {
    if (type == "custom") {
      stop("For `type = \"custom\"`, `model` must be supplied (a named ",
           "list with a `class` field).", call. = FALSE)
    }
    return(default_ng_model(type, fclass))
  }
  if (!is.list(model) || is.null(names(model)) || any(!nzchar(names(model)))) {
    stop("`model` must be NULL or a fully named list.", call. = FALSE)
  }

  if (fclass == "categorical") {
    if (is.null(model$reference)) {
      model$reference <- "marginal"
    }
    if (!is.character(model$reference) || length(model$reference) != 1L ||
        !model$reference %in% c("marginal", "uniform")) {
      stop("`model$reference` must be \"marginal\" or \"uniform\" for a ",
           "categorical feature.", call. = FALSE)
    }
    return(model)
  }

  if (fclass == "continuous") {
    if (is.null(model$reference)) {
      model$reference <- "uniform"
    }
    if (!is.character(model$reference) || length(model$reference) != 1L ||
        !model$reference %in% c("uniform", "empirical")) {
      stop("`model$reference` must be \"uniform\" or \"empirical\" for a ",
           "continuous feature.", call. = FALSE)
    }
    if (model$reference == "uniform") {
      if (is.null(model$range)) {
        model$range <- c(1, 80)
      }
      rg <- model$range
      if (!is.numeric(rg) || length(rg) != 2L || anyNA(rg) ||
          !all(is.finite(rg)) || rg[2L] <= rg[1L]) {
        stop("`model$range` must be a finite increasing numeric pair ",
             "c(lo, hi).", call. = FALSE)
      }
    }
    return(model)
  }

  # fclass == "date"
  if (is.null(model$search)) {
    model$search <- "open"
  }
  if (!is.character(model$search) || length(model$search) != 1L ||
      !model$search %in% c("open", "closed")) {
    stop("`model$search` must be \"open\" or \"closed\" for a date ",
         "feature.", call. = FALSE)
  }
  if (is.null(model$cuts)) {
    model$cuts <- c(-120, -30, 30, 120, 240, 360)
  }
  cuts <- model$cuts
  if (!is.numeric(cuts) || length(cuts) < 1L || anyNA(cuts) ||
      !all(is.finite(cuts)) || is.unsorted(cuts, strictly = TRUE)) {
    stop("`model$cuts` must be a strictly increasing finite numeric ",
         "vector.", call. = FALSE)
  }
  model
}

#' @noRd
validate_ng_freqs <- function(x, tol) {
  if (!is.numeric(x) || is.array(x)) {
    stop("`db_or_freqs` must be a named numeric vector of category ",
         "frequencies for a categorical feature.", call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm)) || anyNA(nm)) {
    stop("`db_or_freqs` must be named (one name per category).",
         call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("`db_or_freqs` must not contain duplicated category labels.",
         call. = FALSE)
  }
  if (length(x) < 2L) {
    stop("`db_or_freqs` must contain at least two categories.",
         call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`db_or_freqs` must not contain NA values.", call. = FALSE)
  }
  if (any(x < 0) || any(x > 1)) {
    stop("`db_or_freqs` entries must lie in [0, 1].", call. = FALSE)
  }
  s <- sum(x)
  if (abs(s - 1) > tol) {
    stop(sprintf("`db_or_freqs` must sum to 1 (within tol=%g); got ",
                 tol), sprintf("sum=%.10f.", s), call. = FALSE)
  }
  x
}

#' @noRd
validate_ng_db <- function(db_or_freqs, fclass, model, tol) {
  if (fclass == "categorical") {
    return(validate_ng_freqs(db_or_freqs, tol))
  }

  if (fclass == "continuous") {
    if (identical(model$reference, "uniform")) {
      if (!is.null(db_or_freqs)) {
        stop("`db_or_freqs` must be NULL when the continuous reference is ",
             "\"uniform\" (the support comes from `model$range`).",
             call. = FALSE)
      }
      return(NULL)
    }
    # empirical
    if (!is.numeric(db_or_freqs) || length(db_or_freqs) < 2L ||
        anyNA(db_or_freqs) || !all(is.finite(db_or_freqs))) {
      stop("`db_or_freqs` must be a finite numeric sample of length >= 2 ",
           "when the continuous reference is \"empirical\".",
           call. = FALSE)
    }
    return(db_or_freqs)
  }

  # fclass == "date"
  nbins <- length(model$cuts) + 1L
  if (identical(model$search, "open")) {
    if (!is.null(db_or_freqs)) {
      stop("`db_or_freqs` must be NULL for an open date search ",
           "(H2 is uniform over the discrepancy bins).", call. = FALSE)
    }
    return(NULL)
  }
  # closed
  if (is.data.frame(db_or_freqs)) {
    if (nrow(db_or_freqs) < 1L) {
      stop("`db_or_freqs` data.frame for a closed date search must have ",
           "at least one row.", call. = FALSE)
    }
    return(db_or_freqs)
  }
  if (!is.numeric(db_or_freqs) || length(db_or_freqs) != nbins ||
      anyNA(db_or_freqs) || any(db_or_freqs < 0) ||
      !all(is.finite(db_or_freqs))) {
    stop(sprintf(paste0("`db_or_freqs` for a closed date search must be a ",
                        "non-negative numeric vector of length %d ",
                        "(length(cuts) + 1) or a data.frame of declared ",
                        "dates."), nbins), call. = FALSE)
  }
  db_or_freqs
}

#' @noRd
validate_ng_error <- function(error, fclass, cats, model, tol) {
  if (fclass == "categorical") {
    if (is.null(error)) {
      error <- 0.05
    }
    if (is.matrix(error)) {
      k <- length(cats)
      if (!is.numeric(error) || nrow(error) != k || ncol(error) != k) {
        stop(sprintf("`error` confusion matrix must be %d x %d ", k, k),
             "(one row/column per category).", call. = FALSE)
      }
      if (anyNA(error) || any(error < 0)) {
        stop("`error` confusion matrix entries must be non-negative and ",
             "non-NA.", call. = FALSE)
      }
      rs <- rowSums(error)
      if (any(abs(rs - 1) > tol)) {
        stop("`error` confusion matrix rows must each sum to 1 (within ",
             "tol).", call. = FALSE)
      }
      # The confusion matrix is positional: row i / column j correspond to
      # the i-th / j-th category in the order of `db_or_freqs`. Dimnames,
      # if present, are informational only (legacy error_matrix_hair() is
      # label-keyed but used positionally by lr_hair_color()).
      return(error)
    }
    if (!is_scalar_finite(error) || error < 0 || error >= 1) {
      stop("`error` must be a finite scalar in [0, 1) or a row-stochastic ",
           "confusion matrix for a categorical feature.", call. = FALSE)
    }
    return(error)
  }

  if (fclass == "continuous") {
    if (is.null(error)) {
      error <- 0.05
    }
    if (!is_scalar_finite(error) || error < 0 || error >= 1) {
      stop("`error` must be a finite scalar in [0, 1) for a continuous ",
           "feature.", call. = FALSE)
    }
    return(error)
  }

  # fclass == "date": Dirichlet alpha vector
  if (is.null(error)) {
    error <- rep(1, length(model$cuts) + 1L)
  }
  nbins <- length(model$cuts) + 1L
  if (!is.numeric(error) || length(error) != nbins || anyNA(error) ||
      !all(is.finite(error)) || any(error <= 0)) {
    stop(sprintf(paste0("`error` must be a strictly positive numeric ",
                        "Dirichlet alpha vector of length %d ",
                        "(length(cuts) + 1)."), nbins), call. = FALSE)
  }
  error
}

#' @noRd
validate_ng_observed <- function(observed, fclass, cats) {
  if (fclass == "categorical") {
    if (length(observed) != 1L || is.na(observed)) {
      stop("`observed` must be a single non-NA category value.",
           call. = FALSE)
    }
    obs_chr <- as.character(observed)
    if (!obs_chr %in% cats) {
      stop("`observed` (\"", obs_chr, "\") is not one of the feature ",
           "categories: ", paste(shQuote(cats), collapse = ", "), ".",
           call. = FALSE)
    }
    return(obs_chr)
  }

  if (fclass == "continuous") {
    if (!is_scalar_finite(observed)) {
      stop("`observed` must be a single finite numeric value for a ",
           "continuous feature.", call. = FALSE)
    }
    return(observed)
  }

  # fclass == "date"
  if (length(observed) != 1L || is.na(observed)) {
    stop("`observed` must be a single non-NA value for a date feature.",
         call. = FALSE)
  }
  if (is.numeric(observed)) {
    if (!is.finite(observed)) {
      stop("`observed` day discrepancy must be finite.", call. = FALSE)
    }
    return(as.numeric(observed))
  }
  parsed <- tryCatch(as.Date(observed),
                      error = function(e) NA,
                      warning = function(w) NA)
  if (length(parsed) != 1L || is.na(parsed)) {
    stop("`observed` must be a Date, a \"YYYY-MM-DD\" string, or a finite ",
         "numeric day discrepancy.", call. = FALSE)
  }
  parsed
}

#' Print method for nongenetic_feature
#'
#' @param x A `nongenetic_feature` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.nongenetic_feature <- function(x, ...) {
  cat("<nongenetic_feature>\n")
  cat("  type          : ", x$type,
      " (", x$feature_class, ")\n", sep = "")
  obs <- if (inherits(x$observed, "Date")) {
    format(x$observed)
  } else {
    as.character(x$observed)
  }
  cat("  observed      : ", obs, "\n", sep = "")

  if (x$feature_class == "categorical") {
    cat("  categories    : ", length(x$categories),
        " (", paste(utils::head(x$categories, 6), collapse = ", "),
        if (length(x$categories) > 6L) ", ..." else "", ")\n", sep = "")
    cat("  reference     : ", x$model$reference, "\n", sep = "")
    if (is.matrix(x$error)) {
      cat("  error         : confusion matrix ",
          nrow(x$error), "x", ncol(x$error), "\n", sep = "")
    } else {
      cat("  error         : eps=", format(x$error, digits = 3), "\n",
          sep = "")
    }
  } else if (x$feature_class == "continuous") {
    cat("  reference     : ", x$model$reference, sep = "")
    if (identical(x$model$reference, "uniform")) {
      cat(" [", x$model$range[1L], ", ", x$model$range[2L], "]", sep = "")
    } else {
      cat(" (n=", length(x$db_or_freqs), ")", sep = "")
    }
    cat("\n")
    cat("  error         : eps=", format(x$error, digits = 3), "\n",
        sep = "")
  } else {
    cat("  search        : ", x$model$search,
        " (", length(x$model$cuts) + 1L, " discrepancy bins)\n", sep = "")
    cat("  error         : Dirichlet alpha (length ",
        length(x$error), ")\n", sep = "")
  }
  invisible(x)
}
