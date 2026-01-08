#' Simulate Likelihood Ratios from Preliminary Investigation Data
#'
#' @description
#' Simulates likelihood ratio (LR) distributions based on non-genetic
#' (preliminary investigation) data such as sex, age, region, height, or
#' birth date. This function generates expected LR distributions under both
#' hypotheses:
#' \itemize{
#'   \item H1 (Related): The unidentified person IS the missing person
#'   \item H2 (Unrelated): The unidentified person is NOT the missing person
#' }
#'
#' @param vartype Character. Type of preliminary investigation variable.
#'   Options: "sex", "region", "age", "height", "birthdate".
#' @param numsims Integer. Number of simulations to perform. Default: 1000.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#' @param int Numeric. Interval parameter for "age" and "height" variables.
#'   Defines the estimation range (e.g., if MP age is 55 and int is 10,
#'   the range is 45-65). Default: 5.
#' @param ErrorRate Numeric (0-1). Error rate for observations. Default: 0.05.
#' @param alphaBdate Numeric vector. Alpha parameters for Dirichlet distribution
#'   used in birthdate LR calculations. Usually frequencies of solved cases
#'   in each category. Default: c(1, 4, 60, 11, 6, 4, 4).
#' @param numReg Integer. Number of regions in the case (for "region" variable).
#'   Default: 6.
#' @param MP Value of the MP's characteristic for closed search. If NULL,
#'   open search is performed. For "sex": "F" or "M"; for "age"/"height":
#'   numeric; for "birthdate": date string; for "region": region ID.
#'   Default: NULL.
#' @param database Data frame. Database of POIs for closed search (when MP
#'   is not NULL). Should have columns matching the variable type (e.g.,
#'   "Sex", "Age", "Height", "Region", "DBD"). Can be output from
#'   \code{\link{sim_poi_prelim}}.
#' @param cuts Numeric vector. Cutoff values for birthdate categories.
#'   Days difference between declared and actual birth dates.
#'   Default: c(-120, -30, 30, 120, 240, 360).
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Unrelated}: LR values simulated under H2 (POI is not MP)
#'     \item \code{Related}: LR values simulated under H1 (POI is MP)
#'   }
#'   Each column contains \code{numsims} values.
#'
#' @details
#' \strong{Open Search (MP = NULL):}
#' Used when it's unknown whether the MP is in the database. LR is computed
#' using general population frequencies as the denominator.
#'
#' \strong{Closed Search (MP specified):}
#' Used when comparing a specific MP against a database. LR denominator
#' uses frequencies from the actual database.
#'
#' \strong{Variable-specific calculations:}
#' \itemize{
#'   \item \emph{sex}: Binary match/mismatch with error rate
#'   \item \emph{region}: Match against numReg possible regions
#'   \item \emph{age/height}: Match if within +/- int of MP value
#'   \item \emph{birthdate}: Dirichlet-based probability for date discrepancy
#' }
#'
#' @seealso
#' \code{\link{sim_lr_genetic}} for genetic LR simulations,
#' \code{\link{lr_combine}} for combining LRs from different sources,
#' \code{\link{sim_poi_prelim}} for creating preliminary databases.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Open search for sex variable
#' lr_sex <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#' head(lr_sex)
#'
#' # Check distribution
#' summary(log10(lr_sex$Related))
#' summary(log10(lr_sex$Unrelated))
#'
#' # Visualize
#' plot_lr_distribution(lr_sex)
#'
#' # Closed search with database
#' db <- sim_poi_prelim(numsims = 100, seed = 456)
#' lr_sex_closed <- sim_lr_prelim(
#'   "sex",
#'   numsims = 500,
#'   MP = "F",
#'   database = db
#' )
#'
#' # Age variable
#' lr_age <- sim_lr_prelim("age", numsims = 500, int = 10)

sim_lr_prelim <- function(vartype,
                          numsims = 1000,
                          seed = 123,
                          int = 5,
                          ErrorRate = 0.05,
                          alphaBdate = c(1, 4, 60, 11, 6, 4, 4),
                          numReg = 6,
                          MP = NULL,
                          database = NULL,
                          cuts = c(-120, -30, 30, 120, 240, 360)) {

  # Input validation
  allowed_vartypes <- c("sex", "region", "age", "height", "birthdate")
  if (!vartype %in% allowed_vartypes) {
    stop("vartype must be one of: ", paste(allowed_vartypes, collapse = ", "))
  }
  if (!is.numeric(numsims) || numsims < 1) {
    stop("numsims must be a positive integer")
  }
  if (!is.numeric(ErrorRate) || ErrorRate < 0 || ErrorRate > 1) {
    stop("ErrorRate must be between 0 and 1")
  }
  if (!is.numeric(numReg) || numReg < 2) {
    stop("numReg must be at least 2")
  }

  # Helper function: Dirichlet probability estimation via method of moments
  # Returns list of estimated category probabilities
  .estimate_dirichlet_probs <- function(alpha, n_samples = 500) {
    # Generate Dirichlet samples via gamma distribution
    x <- matrix(stats::rgamma(n_samples * length(alpha), alpha, 1),
                ncol = length(alpha), byrow = TRUE)
    x <- x / rowSums(x)

    # Add small constant to avoid log(0)
    x <- pmax(x, 1e-10)

    # Method of moments estimation
    x1_mean <- mean(x[, 1])
    x1_mean_sq <- mean(x[, 1]^2)
    denom <- x1_mean_sq - x1_mean^2

    # Guard against division by zero in moment estimation
    if (abs(denom) < 1e-10) {
      warning("Dirichlet moment estimation: near-zero denominator, using column means")
      mom <- colMeans(x)
    } else {
      mom <- colMeans(x) * (x1_mean - x1_mean_sq) / denom
    }

    as.list(mom / sum(mom))
  }

  # Helper function: Safe LR calculation with division by zero protection
  .safe_lr <- function(numerator, denominator, varname = "variable") {
    if (abs(denominator) < 1e-10) {
      warning(sprintf("Division by near-zero in %s LR calculation (denom=%.6f). ",
                      varname, denominator),
              "Consider different parameters or larger database.")
      return(ifelse(numerator > 0, 1e6, 1))  # Cap at 1e6 to avoid Inf
    }
    numerator / denominator
  }

  set.seed(seed)

  # Initialize result vectors
  a <- NULL  # Related (H1) LRs
  b <- NULL  # Unrelated (H2) LRs

  if (is.null(MP)) {
    # === OPEN SEARCH ===

    if (vartype == "sex") {
      sexLRvalues <- c((1 - ErrorRate) / 0.5, ErrorRate / 0.5)
      a <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))
    }
    else if (vartype == "region") {
      # P(match | H2) = 1/numReg (uniform probability)
      # P(mismatch | H2) = (numReg-1)/numReg
      p_match_H2 <- 1 / numReg
      p_mismatch_H2 <- (numReg - 1) / numReg
      # LR(match) = P(match|H1) / P(match|H2) = (1-ErrorRate) / (1/numReg)
      # LR(mismatch) = P(mismatch|H1) / P(mismatch|H2) = ErrorRate / ((numReg-1)/numReg)
      RegLRvalues <- c((1 - ErrorRate) / p_match_H2, ErrorRate / p_mismatch_H2)
      a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(p_match_H2, p_mismatch_H2))
    }
    else if (vartype == "age") {
      AgeLRvalues <- c((1 - ErrorRate) / 0.5, ErrorRate / 0.5)
      a <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))
    }
    else if (vartype == "height") {
      heightLRvalues <- c((1 - ErrorRate) / 0.5, ErrorRate / 0.5)
      a <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))
    }
    else if (vartype == "birthdate") {
      # Estimate category probabilities using Dirichlet model
      fit2 <- .estimate_dirichlet_probs(alphaBdate)
      n_cats <- length(fit2)

      # LR values for each category: P(cat|H1) / P(cat|H2_uniform)
      # Under H2, assume uniform distribution across categories
      h2_uniform <- 1 / n_cats
      bdateLRvalues <- lapply(fit2, function(p) p / h2_uniform)

      a <- sample(unlist(bdateLRvalues), numsims, replace = TRUE, prob = unlist(fit2))
      b <- sample(unlist(bdateLRvalues), numsims, replace = TRUE)
    }
  }
  else {
    # === CLOSED SEARCH ===

    if (is.null(database)) {
      stop("database parameter is required for closed search (when MP is specified)")
    }

    if (vartype == "sex") {
      if (!"Sex" %in% names(database)) {
        stop("database must have 'Sex' column for sex variable")
      }
      fs <- sum(database$Sex == MP) / length(database$Sex)
      # Safe LR calculation with division by zero protection
      lr_match <- .safe_lr(1 - ErrorRate, fs, "sex match")
      lr_mismatch <- .safe_lr(ErrorRate, 1 - fs, "sex mismatch")
      sexLRvalues <- c(lr_match, lr_mismatch)
      # Ensure valid probabilities for sampling
      fs_safe <- max(min(fs, 1 - 1e-10), 1e-10)
      a <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(fs_safe, 1 - fs_safe))
    }
    else if (vartype == "region") {
      if (!"Region" %in% names(database)) {
        stop("database must have 'Region' column for region variable")
      }
      fr <- sum(database$Region == MP) / length(database$Region)
      lr_match <- .safe_lr(1 - ErrorRate, fr, "region match")
      lr_mismatch <- .safe_lr(ErrorRate, 1 - fr, "region mismatch")
      RegLRvalues <- c(lr_match, lr_mismatch)
      fr_safe <- max(min(fr, 1 - 1e-10), 1e-10)
      a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(fr_safe, 1 - fr_safe))
    }
    else if (vartype == "age") {
      if (!"Age" %in% names(database)) {
        stop("database must have 'Age' column for age variable")
      }
      fa <- sum(database$Age < (MP + int) & database$Age > (MP - int)) / length(database$Age)
      lr_match <- .safe_lr(1 - ErrorRate, fa, "age match")
      lr_mismatch <- .safe_lr(ErrorRate, 1 - fa, "age mismatch")
      AgeLRvalues <- c(lr_match, lr_mismatch)
      fa_safe <- max(min(fa, 1 - 1e-10), 1e-10)
      a <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(fa_safe, 1 - fa_safe))
    }
    else if (vartype == "height") {
      if (!"Height" %in% names(database)) {
        stop("database must have 'Height' column for height variable")
      }
      fh <- sum(database$Height < (MP + int) & database$Height > (MP - int)) / length(database$Height)
      lr_match <- .safe_lr(1 - ErrorRate, fh, "height match")
      lr_mismatch <- .safe_lr(ErrorRate, 1 - fh, "height mismatch")
      heightLRvalues <- c(lr_match, lr_mismatch)
      fh_safe <- max(min(fh, 1 - 1e-10), 1e-10)
      a <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(fh_safe, 1 - fh_safe))
    }
    else if (vartype == "birthdate") {
      if (!"DBD" %in% names(database)) {
        stop("database must have 'DBD' column for birthdate variable")
      }
      MP <- as.Date(MP)
      PrelimData <- data.frame(
        ABD = MP,
        DBD = database$DBD,
        Dis = as.numeric(database$DBD - MP)
      )

      cut_results <- cut(PrelimData$Dis, breaks = c(-Inf, cuts, Inf))
      freq_table <- table(cut_results)
      relfreq <- as.vector(freq_table / sum(freq_table))

      # Add pseudocount to avoid zero frequencies
      alpha2 <- pmax(relfreq, 1e-6)
      alpha2 <- alpha2 / sum(alpha2)  # Renormalize

      # Estimate H1 probabilities using Dirichlet model
      fit2 <- .estimate_dirichlet_probs(alphaBdate)

      # LR values using database frequencies with safe division
      bdateLRvalues <- lapply(seq_along(fit2), function(i) {
        .safe_lr(fit2[[i]], alpha2[i], paste0("birthdate cat ", i))
      })

      a <- sample(unlist(bdateLRvalues), numsims, replace = TRUE, prob = unlist(fit2))
      b <- sample(unlist(bdateLRvalues), numsims, replace = TRUE, prob = alpha2)
    }
  }

  # Return as dataframe
  LRsimulated <- cbind(b, a)
  colnames(LRsimulated) <- c("Unrelated", "Related")
  return(as.data.frame(LRsimulated))
}
