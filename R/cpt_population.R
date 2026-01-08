#' Population-Based Conditional Probability Table
#'
#' @description
#' Computes a conditional probability table (CPT) representing the joint
#' probability distribution of sex, age group, and hair color in the reference
#' population. This table represents P(D|H2), the probability of observing
#' the evidence under the hypothesis that the unidentified person is NOT
#' the missing person.
#'
#' The function assumes a uniform age distribution across the population
#' and computes the probability of falling within or outside the specified
#' age range.
#'
#' @param propS Numeric vector of length 2. Sex proportions in the population,
#'   specified as c(proportion_female, proportion_male). Must sum to 1.
#'   Default: c(0.5, 0.5) for equal proportions.
#' @param MPa Numeric. Missing person's estimated age in years. Used to define
#'   the center of the age matching interval. Default: 40.
#' @param MPr Numeric. Age range tolerance in years (plus/minus). The matching
#'   age interval is \code{(MPa - MPr)} to \code{(MPa + MPr)}. Individuals within
#'   this range are classified as T1 (age match), others as T0 (age mismatch).
#'   Default: 6.
#' @param propC Numeric vector of length 5. Hair color proportions in the
#'   population for colors 1 through 5. Must sum to 1. Default values represent
#'   a typical distribution. Colors are: 1=Black, 2=Brown, 3=Blonde,
#'   4=Red, 5=Gray/White.
#'
#' @return A 4x5 numeric matrix representing joint probabilities. Rows
#'   correspond to sex-age group combinations:
#'   \itemize{
#'     \item F-T1: Female, age within range
#'     \item F-T0: Female, age outside range
#'     \item M-T1: Male, age within range
#'     \item M-T0: Male, age outside range
#'   }
#'   Columns correspond to hair colors 1-5. Each cell contains
#'   P(Sex, AgeGroup, HairColor | H2).
#'
#' @details
#' The probability of age match (T1) is calculated assuming a uniform
#' distribution over ages 1-80:
#' \deqn{P(T1) = (MPa + MPr - (MPa - MPr)) / 80 = 2 \times MPr / 80}
#'
#' The joint probability for each cell is:
#' \deqn{P(Sex, Age, Color) = P(Sex) \times P(AgeGroup) \times P(Color)}
#'
#' @seealso
#' \code{\link{cpt_missing_person}} for the H1 conditional probability table,
#' \code{\link{plot_cpt}} for visualization of CPTs.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Default parameters: equal sex proportions, MP age 40 +/- 6 years
#' cpt_h2 <- cpt_population()
#' print(cpt_h2)
#'
#' # Custom population: 60% female, narrower age range, different hair colors
#' cpt_custom <- cpt_population(
#'   propS = c(0.6, 0.4),
#'   MPa = 35,
#'   MPr = 3,
#'   propC = c(0.4, 0.25, 0.2, 0.1, 0.05)
#' )
#'
#' # Verify rows sum correctly (should sum to hair color proportions)
#' colSums(cpt_h2)

cpt_population <- function(propS = c(0.5, 0.5),
                           MPa = 40,
                           MPr = 6,
                           propC = c(0.3, 0.2, 0.25, 0.15, 0.1)) {

  # Input validation
  if (!is.numeric(propS) || length(propS) != 2) {
    stop("propS must be a numeric vector of length 2")
  }
  if (any(propS < 0) || any(propS > 1)) {
    stop("propS values must be between 0 and 1")
  }
  if (abs(sum(propS) - 1) > 1e-6) {
    warning("propS does not sum to 1; normalizing")
    propS <- propS / sum(propS)
  }

  if (!is.numeric(MPa) || length(MPa) != 1 || MPa < 1 || MPa > 80) {
    stop("MPa must be a numeric value between 1 and 80")
  }
  if (!is.numeric(MPr) || length(MPr) != 1 || MPr < 0) {
    stop("MPr must be a non-negative numeric value")
  }
  if (MPa - MPr < 1 || MPa + MPr > 80) {
    warning("Age range extends beyond [1, 80]; clamping to valid range")
  }

  if (!is.numeric(propC) || length(propC) != 5) {
    stop("propC must be a numeric vector of length 5")
  }
  if (any(propC < 0)) {
    stop("propC values must be non-negative")
  }
  if (abs(sum(propC) - 1) > 1e-6) {
    warning("propC does not sum to 1; normalizing")
    propC <- propC / sum(propC)
  }

  # Age range for population (assumed uniform 1-80)
  Age <- seq(1, 80)

  # Calculate age boundaries

  MPmin <- MPa - MPr
  MPmax <- MPa + MPr

  # Probability of falling within age range (T1) assuming uniform distribution

  T1p <- (MPmax - MPmin) / length(Age)
  T0p <- 1 - T1p

  # Create joint probability vector for sex-age combinations
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
  jointprob <- c(
    propS[1] * T1p,  # Female, age match
    propS[1] * T0p,  # Female, age mismatch
    propS[2] * T1p,  # Male, age match
    propS[2] * T0p   # Male, age mismatch
  )
  names(jointprob) <- jointname

  # Create CPT as outer product with hair color proportions
  CPTpop <- outer(jointprob, propC)
  colnames(CPTpop) <- as.character(1:5)


  return(CPTpop)
}
