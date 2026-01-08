#' Likelihood Ratio for Birth Date
#'
#' @description
#' Computes likelihood ratios (LRs) based on the discrepancy between the
#' actual birth date (ABD) of the missing person and the declared birth
#' date (DBD) of the person of interest. Uses Dirichlet distribution to
#' model category probabilities.
#'
#' @param ABD Character or Date. Actual birth date of the missing person
#'   in "YYYY-MM-DD" format. Default: "1976-05-31".
#' @param DBD Character or Date. Declared birth date of the person of
#'   interest in "YYYY-MM-DD" format. Default: "1976-07-15".
#' @param alpha Numeric vector. Alpha parameters for the Dirichlet distribution,
#'   typically representing frequencies of solved cases in each discrepancy
#'   category. Length should be one more than length of \code{cuts}.
#'   Default: c(1, 4, 60, 11, 6, 4, 4).
#' @param cuts Numeric vector. Cutoff values (in days) for categorizing the
#'   difference between DBD and ABD. Creates length(cuts)+1 categories.
#'   Default: c(-120, -30, 30, 120, 240, 360).
#' @param type Integer (1 or 2). Type of search scenario:
#'   \itemize{
#'     \item 1: Open search - MP may not be in database (uses uniform H2)
#'     \item 2: Closed search - MP is in database (uses database frequencies)
#'   }
#'   Default: 1.
#' @param PrelimData Data.frame. Required when type = 2. Contains DBD column
#'   for persons of interest in the database. Can be output from
#'   \code{\link{sim_poi_prelim}}.
#' @param draw Integer. Number of Dirichlet samples for probability estimation.
#'   Default: 500.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#'
#' @return Numeric. The likelihood ratio for the given birth date discrepancy.
#'   Also printed to console.
#'
#' @details
#' \strong{Categories:}
#' The difference between DBD and ABD (in days) is categorized using the
#' \code{cuts} vector. Default categories are:
#' \enumerate{
#'   \item < -120 days (DBD more than 4 months before ABD)
#'   \item -120 to -30 days
#'   \item -30 to 30 days (close match)
#'   \item 30 to 120 days
#'   \item 120 to 240 days
#'   \item 240 to 360 days
#'   \item > 360 days (DBD more than 1 year after ABD)
#' }
#'
#' \strong{Dirichlet Model:}
#' Uses method of moments to estimate category probabilities from Dirichlet
#' samples. The \code{alpha} parameter reflects prior knowledge from solved
#' cases about the distribution of birth date discrepancies.
#'
#' \strong{LR Calculation:}
#' \itemize{
#'   \item Type 1: LR = P(category | H1) / (1/n_categories)
#'   \item Type 2: LR = P(category | H1) / P(category in database)
#' }
#'
#' @seealso
#' \code{\link{sim_lr_prelim}} for simulating LR distributions,
#' \code{\link{sim_poi_prelim}} for generating preliminary databases.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @import DirichletReg
#' @import dplyr
#' @examples
#' # Type 1: Open search - close match (45 days difference)
#' lr1 <- lr_birthdate(
#'   ABD = "1976-05-31",
#'   DBD = "1976-07-15",
#'   type = 1,
#'   seed = 123
#' )
#'
#' # Type 1: Open search - larger discrepancy
#' lr2 <- lr_birthdate(
#'   ABD = "1976-05-31",
#'   DBD = "1977-03-15",
#'   type = 1,
#'   seed = 123
#' )
#'
#' \dontrun{
#' # Type 2: Closed search with database
#' # Requires a large database with varied birth dates
#' db <- sim_poi_prelim(numsims = 1000, seed = 456)
#' lr3 <- lr_birthdate(
#'   ABD = "1976-05-31",
#'   DBD = "1976-07-15",
#'   type = 2,
#'   PrelimData = db,
#'   seed = 123
#' )
#' }

lr_birthdate <- function(ABD = "1976-05-31",
                         DBD = "1976-07-15",
                         alpha = c(1, 4, 60, 11, 6, 4, 4),
                         cuts = c(-120, -30, 30, 120, 240, 360),
                         type = 1,
                         PrelimData = NULL,
                         draw = 500,
                         seed = 123) {

  # Input validation
  if (!type %in% c(1, 2)) {
    stop("type must be 1 (open search) or 2 (closed search)")
  }

  if (length(alpha) != length(cuts) + 1) {
    stop("length(alpha) must equal length(cuts) + 1. ",
         "alpha has ", length(alpha), " elements, cuts has ", length(cuts),
         " elements (expected ", length(cuts) + 1, " alpha values).")
  }

  if (any(alpha <= 0)) {
    stop("All alpha values must be positive")
  }

  if (draw < 10) {
    stop("draw must be at least 10 for reliable estimation")
  }

  # Validate date formats
  ABD <- tryCatch(as.Date(ABD), error = function(e) {
    stop("ABD must be a valid date in 'YYYY-MM-DD' format")
  })
  DBD <- tryCatch(as.Date(DBD), error = function(e) {
    stop("DBD must be a valid date in 'YYYY-MM-DD' format")
  })

  set.seed(seed)

  # Generate Dirichlet samples for H1 probabilities
  x <- DirichletReg::rdirichlet(draw, alpha)

  # Method of moments estimation for Dirichlet parameters
  # Uses the expected values of the Dirichlet samples
  mom <- colMeans(x) * (mean(x[, 1]) - mean(x[, 1]^2)) /
         (mean(x[, 1]^2) - (mean(x[, 1]))^2)
  fit2 <- as.list(mom / sum(mom))

  # Calculate difference in days
  Dis0 <- julian(DBD, ABD)
  Dist <- Dis0[1]

  # Determine which category the discrepancy falls into
  # Categories: 1 = below cuts[1], 2 = cuts[1] to cuts[2], ..., n = above cuts[n-1]
  n_cats <- length(alpha)
  n_cuts <- length(cuts)

  if (Dist < cuts[1]) {
    H1 <- fit2[[1]]
  } else if (Dist >= cuts[n_cuts]) {
    H1 <- fit2[[n_cats]]
  } else {
    # Find which interval Dist falls into
    for (i in seq_len(n_cuts - 1)) {
      if (Dist >= cuts[i] && Dist < cuts[i + 1]) {
        H1 <- fit2[[i + 1]]
        break
      }
    }
  }

  if (type == 1) {
    # Open search: uniform H2
    H2 <- 1 / length(alpha)
    LR <- as.numeric(H1) / H2
  }

  if (type == 2) {
    # Closed search: use database frequencies
    if (is.null(PrelimData)) {
      stop("PrelimData is required for type = 2 (closed search)")
    }

    PrelimData <- dplyr::mutate(PrelimData, Dis = julian(DBD, ABD))
    PrelimData <- as.data.frame(table(cut(PrelimData$Dis, breaks = c(-Inf, cuts, Inf))))
    alpha2 <- as.vector(PrelimData$Freq)

    # Generate Dirichlet samples for H2
    x2 <- DirichletReg::rdirichlet(draw, alpha2)
    temp2 <- dim(x2)
    n <- temp2[1]
    m2 <- temp2[2]
    lpb2 <- apply(log(x2), 2, mean)
    mom2 <- apply(x2, 2, mean) * (mean(x2[, 1]) - mean(x2[, 1]^2)) / (mean(x2[, 1]^2) - ((mean(x2[, 1]))^2))
    fit3 <- as.list(mom2 / sum(mom2))

    # Find H2 probability for the category
    n_cats2 <- length(alpha2)
    if (Dist < cuts[1]) {
      H2 <- fit3[[1]]
    } else if (Dist >= cuts[n_cuts]) {
      H2 <- fit3[[n_cats2]]
    } else {
      for (i in seq_len(n_cuts - 1)) {
        if (Dist >= cuts[i] && Dist < cuts[i + 1]) {
          H2 <- fit3[[i + 1]]
          break
        }
      }
    }

    LR <- as.numeric(H1) / as.numeric(H2)
  }

  message(paste("LR =", round(LR, 4)))
  invisible(LR)
}
