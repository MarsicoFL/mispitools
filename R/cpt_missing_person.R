#' Missing Person-Based Conditional Probability Table
#'
#' @description
#' Computes a conditional probability table (CPT) representing the probability
#' of observing evidence given the hypothesis that the unidentified person IS
#' the missing person. This table represents P(D|H1), accounting for potential
#' observation errors in sex, age, and hair color.
#'
#' The function incorporates error rates (epsilon values) that model the
#' probability of misclassifying the true characteristics of the missing
#' person during observation.
#'
#' @param MPs Character. Missing person's biological sex: "F" for female,
#'   "M" for male. Default: "F".
#' @param MPc Integer (1-5). Missing person's hair color category:
#'   1=Black, 2=Brown, 3=Blonde, 4=Red, 5=Gray/White. Default: 1.
#' @param eps Numeric (0-1). Error rate for sex observation. The probability
#'   of incorrectly recording the sex. Default: 0.05.
#' @param epa Numeric (0-1). Error rate for age categorization. The probability
#'   of classifying a person in the wrong age group (T0 instead of T1).
#'   Default: 0.05.
#' @param epc Matrix. Hair color error/confusion matrix, typically created
#'   with \code{\link{error_matrix_hair}}. Rows represent true colors,
#'   columns represent observed colors. Default: \code{error_matrix_hair()}.
#'
#' @return A 4x5 numeric matrix representing conditional probabilities under H1.
#'   Rows correspond to observed sex-age group combinations:
#'   \itemize{
#'     \item F-T1: Observed as Female, age within range
#'     \item F-T0: Observed as Female, age outside range
#'     \item M-T1: Observed as Male, age within range
#'     \item M-T0: Observed as Male, age outside range
#'   }
#'   Columns correspond to observed hair colors 1-5.
#'   Each cell contains P(Observed Sex, Observed Age, Observed Color | H1, MP characteristics).
#'
#' @details
#' For a female MP (MPs = "F"), the joint probabilities are:
#' \itemize{
#'   \item P(F-T1) = (1 - eps) * (1 - epa): Correctly observed sex and age
#'   \item P(F-T0) = (1 - eps) * epa: Correct sex, wrong age group
#'   \item P(M-T1) = eps * (1 - epa): Wrong sex, correct age
#'   \item P(M-T0) = eps * epa: Wrong sex and age
#' }
#'
#' The hair color probabilities come from the error matrix row corresponding
#' to the MP's true hair color.
#'
#' @seealso
#' \code{\link{cpt_population}} for the H2 conditional probability table,
#' \code{\link{error_matrix_hair}} for creating the color error matrix,
#' \code{\link{plot_cpt}} for visualization of CPTs.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Default: Female MP with black hair
#' cpt_h1 <- cpt_missing_person()
#' print(cpt_h1)
#'
#' # Male MP with brown hair, higher error rates
#' cpt_h1_male <- cpt_missing_person(
#'   MPs = "M",
#'   MPc = 2,
#'   eps = 0.10,
#'   epa = 0.10
#' )
#'
#' # Compare H1 and H2 to compute LR
#' cpt_h2 <- cpt_population()
#' lr_matrix <- cpt_h1 / cpt_h2
#' print(log10(lr_matrix))

cpt_missing_person <- function(MPs = "F",
                                MPc = 1,
                                eps = 0.05,
                                epa = 0.05,
                                epc = error_matrix_hair()) {

  # Handle numeric input for sex (1=Female, 2=Male)
  if (is.numeric(MPs)) {
    if (!MPs %in% c(1, 2)) {
      stop("Numeric MPs must be 1 (Female) or 2 (Male)")
    }
    MPs <- if (MPs == 1) "F" else "M"
  }

  # Input validation
  if (!MPs %in% c("F", "M")) {
    stop("MPs must be 'F' (Female) or 'M' (Male)")
  }

  if (!is.numeric(MPc) || length(MPc) != 1 || !MPc %in% 1:5) {
    stop("MPc must be an integer from 1 to 5 (hair color category)")
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps < 0 || eps > 1) {
    stop("eps must be a numeric value between 0 and 1")
  }

  if (!is.numeric(epa) || length(epa) != 1 || epa < 0 || epa > 1) {
    stop("epa must be a numeric value between 0 and 1")
  }

  if (!is.matrix(epc) || nrow(epc) != 5 || ncol(epc) != 5) {
    stop("epc must be a 5x5 matrix (hair color error matrix)")
  }

  # Check that error matrix rows sum to approximately 1
  row_sums <- rowSums(epc)
  if (any(abs(row_sums - 1) > 1e-6)) {
    warning("Error matrix rows do not sum to 1; normalizing")
    epc <- epc / row_sums
  }

  # Create joint probability vector for sex-age combinations
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")

  if (MPs == "F") {
    # Female MP: high prob of observing female, low prob of male
    jointprob <- c(
      (1 - eps) * (1 - epa),  # F-T1: correct sex, correct age group
      (1 - eps) * epa,        # F-T0: correct sex, wrong age group
      eps * (1 - epa),        # M-T1: wrong sex, correct age
      eps * epa               # M-T0: wrong sex, wrong age
    )
  } else if (MPs == "M") {
    # Male MP: high prob of observing male, low prob of female
    jointprob <- c(
      eps * (1 - epa),        # F-T1: wrong sex, correct age
      eps * epa,              # F-T0: wrong sex, wrong age
      (1 - eps) * (1 - epa),  # M-T1: correct sex, correct age group
      (1 - eps) * epa         # M-T0: correct sex, wrong age group
    )
  }
  names(jointprob) <- jointname

  # Get hair color probabilities from error matrix (row = true color)
  Col <- c(1, 2, 3, 4, 5)
  probC <- epc[MPc, ]
  names(probC) <- Col

  # Create CPT as outer product
  CPTmp <- outer(jointprob, probC)
  colnames(CPTmp) <- as.character(1:5)

  return(CPTmp)
}
