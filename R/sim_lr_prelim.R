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
      H2 <- 1 / numReg
      RegLRvalues <- c((1 - ErrorRate) / (1 / H2), ErrorRate / H2)
      a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(RegLRvalues, numsims, replace = TRUE)
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
      # Generate Dirichlet samples
      x <- matrix(stats::rgamma(500 * length(alphaBdate), alphaBdate, 1),
                  ncol = length(alphaBdate), byrow = TRUE)
      x <- x / rowSums(x)

      # Method of moments estimation
      lpb <- colMeans(log(x))
      mom <- colMeans(x) * (mean(x[, 1]) - mean(x[, 1]^2)) / (mean(x[, 1]^2) - (mean(x[, 1]))^2)
      fit2 <- as.list(mom / sum(mom))

      # LR values for each category
      bdateLRvalues <- lapply(seq_along(fit2), function(i) fit2[[i]] / (1 / length(fit2)))

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
      fs <- sum(database$Sex == MP) / length(database$Sex)
      sexLRvalues <- c((1 - ErrorRate) / fs, ErrorRate / (1 - fs))
      a <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(fs, 1 - fs))
    }
    else if (vartype == "region") {
      fr <- sum(database$Region == MP) / length(database$Region)
      RegLRvalues <- c((1 - ErrorRate) / fr, ErrorRate / (1 - fr))
      a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(fr, 1 - fr))
    }
    else if (vartype == "age") {
      fa <- sum(database$Age < (MP + int) & database$Age > (MP - int)) / length(database$Age)
      AgeLRvalues <- c((1 - ErrorRate) / fa, ErrorRate / (1 - fa))
      a <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(fa, 1 - fa))
    }
    else if (vartype == "height") {
      fh <- sum(database$Height < (MP + int) & database$Height > (MP - int)) / length(database$Height)
      heightLRvalues <- c((1 - ErrorRate) / fh, ErrorRate / (1 - fh))
      a <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(1 - ErrorRate, ErrorRate))
      b <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(fh, 1 - fh))
    }
    else if (vartype == "birthdate") {
      MP <- as.Date(MP)
      PrelimData <- data.frame(
        ABD = MP,
        DBD = database$DBD,
        Dis = as.numeric(database$DBD - MP)
      )

      cut_results <- cut(PrelimData$Dis, breaks = c(-Inf, cuts, Inf))
      freq_table <- table(cut_results)
      relfreq <- as.vector(freq_table / sum(freq_table))
      alpha2 <- relfreq

      # Generate Dirichlet samples
      x <- matrix(stats::rgamma(500 * length(alphaBdate), alphaBdate, 1),
                  ncol = length(alphaBdate), byrow = TRUE)
      x <- x / rowSums(x)

      # Method of moments estimation
      lpb <- colMeans(log(x))
      mom <- colMeans(x) * (mean(x[, 1]) - mean(x[, 1]^2)) / (mean(x[, 1]^2) - (mean(x[, 1]))^2)
      fit2 <- as.list(mom / sum(mom))

      # LR values using database frequencies
      bdateLRvalues <- lapply(seq_along(fit2), function(i) fit2[[i]] / alpha2[i])

      a <- sample(unlist(bdateLRvalues), numsims, replace = TRUE, prob = unlist(fit2))
      b <- sample(unlist(bdateLRvalues), numsims, replace = TRUE, prob = alpha2)
    }
  }

  # Return as dataframe
  LRsimulated <- cbind(b, a)
  colnames(LRsimulated) <- c("Unrelated", "Related")
  return(as.data.frame(LRsimulated))
}
