#' Hair Color Error/Confusion Matrix
#'
#' @description
#' Creates a 5x5 error matrix (also known as confusion matrix) that models
#' the probability of observing each hair color given the true hair color.
#' This accounts for observation errors in hair color classification.
#'
#' The matrix rows represent the true hair color of the missing person,
#' and columns represent the observed hair color. Each row sums to 1,
#' indicating that some color must be observed.
#'
#' @param errorModel Character. Type of error model to use:
#'   \itemize{
#'     \item "custom": Use specific error rates for each color pair (default)
#'     \item "uniform": Use a single error rate for all color pairs
#'   }
#' @param ep Numeric (0-1). Base error rate used when errorModel = "uniform".
#'   Represents the probability of confusing any two different colors.
#'   Default: 0.01.
#' @param ep12 Numeric. Error rate between colors 1 (Black) and 2 (Brown).
#'   Default: 0.01.
#' @param ep13 Numeric. Error rate between colors 1 (Black) and 3 (Blonde).
#'   Default: 0.005.
#' @param ep14 Numeric. Error rate between colors 1 (Black) and 4 (Red).
#'   Default: 0.01.
#' @param ep15 Numeric. Error rate between colors 1 (Black) and 5 (Gray/White).
#'   Default: 0.003.
#' @param ep23 Numeric. Error rate between colors 2 (Brown) and 3 (Blonde).
#'   Default: 0.01.
#' @param ep24 Numeric. Error rate between colors 2 (Brown) and 4 (Red).
#'   Default: 0.003.
#' @param ep25 Numeric. Error rate between colors 2 (Brown) and 5 (Gray/White).
#'   Default: 0.01.
#' @param ep34 Numeric. Error rate between colors 3 (Blonde) and 4 (Red).
#'   Default: 0.003.
#' @param ep35 Numeric. Error rate between colors 3 (Blonde) and 5 (Gray/White).
#'   Default: 0.003.
#' @param ep45 Numeric. Error rate between colors 4 (Red) and 5 (Gray/White).
#'   Default: 0.01.
#'
#' @return A 5x5 numeric matrix where:
#'   \itemize{
#'     \item Rows represent true hair colors (1-5)
#'     \item Columns represent observed hair colors (1-5)
#'     \item Cell (i,j) = P(observed color j | true color i)
#'     \item Each row sums to 1
#'     \item Diagonal elements are highest (correct observations)
#'   }
#'   Hair color codes: 1=Black, 2=Brown, 3=Blonde, 4=Red, 5=Gray/White.
#'
#' @details
#' The error rates are symmetric: the probability of confusing color A with
#' color B equals the probability of confusing B with A.
#'
#' The diagonal elements (correct observations) are calculated to ensure
#' each row sums to 1:
#' \deqn{P(i|i) = 1 / (1 + \sum_{j \neq i} ep_{ij})}
#'
#' Lower error rates between dissimilar colors (e.g., black and blonde)
#' and higher rates between similar colors (e.g., brown and red) reflect
#' realistic observation patterns.
#'
#' @seealso
#' \code{\link{cpt_missing_person}} which uses this matrix,
#' \code{\link{lr_hair_color}} for hair color LR calculations.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Default custom error model
#' emat <- error_matrix_hair()
#' print(round(emat, 4))
#'
#' # Verify rows sum to 1
#' rowSums(emat)
#'
#' # Uniform error model with 2% error rate
#' emat_uniform <- error_matrix_hair(errorModel = "uniform", ep = 0.02)
#' print(round(emat_uniform, 4))
#'
#' # Higher error rates for similar colors
#' emat_custom <- error_matrix_hair(
#'   errorModel = "custom",
#'   ep12 = 0.05,  # Black-Brown confusion more likely
#'   ep23 = 0.05,  # Brown-Blonde confusion more likely
#'   ep34 = 0.05   # Blonde-Red confusion more likely
#' )

error_matrix_hair <- function(errorModel = c("custom", "uniform")[1],
                               ep = 0.01,
                               ep12 = 0.01,
                               ep13 = 0.005,
                               ep14 = 0.01,
                               ep15 = 0.003,
                               ep23 = 0.01,
                               ep24 = 0.003,
                               ep25 = 0.01,
                               ep34 = 0.003,
                               ep35 = 0.003,
                               ep45 = 0.01) {

  # If uniform model, use same error rate for all pairs
  if (errorModel == "uniform") {
    ep12 <- ep13 <- ep14 <- ep15 <- ep
    ep23 <- ep24 <- ep25 <- ep
    ep34 <- ep35 <- ep
    ep45 <- ep
  }

  # Calculate diagonal elements (normalization factors)
  # Each ensures the row sums to 1
  l1 <- 1 / (1 + ep12 + ep13 + ep14 + ep15)
  l2 <- 1 / (1 + ep12 + ep23 + ep24 + ep25)
  l3 <- 1 / (1 + ep13 + ep23 + ep34 + ep35)
  l4 <- 1 / (1 + ep14 + ep24 + ep34 + ep45)
  l5 <- 1 / (1 + ep15 + ep25 + ep35 + ep45)

  # Build error matrix
  # Row i = probabilities of observing each color given true color i
  errorMat <- rbind(
    c(l1,       l1*ep12,  l1*ep13,  l1*ep14,  l1*ep15),  # True color 1 (Black)
    c(l2*ep12,  l2,       l2*ep23,  l2*ep24,  l2*ep25),  # True color 2 (Brown)
    c(l3*ep13,  l3*ep23,  l3,       l3*ep34,  l3*ep35),  # True color 3 (Blonde)
    c(l4*ep14,  l4*ep24,  l4*ep34,  l4,       l4*ep45),  # True color 4 (Red)
    c(l5*ep15,  l5*ep25,  l5*ep35,  l5*ep45,  l5)        # True color 5 (Gray/White)
  )

  rownames(errorMat) <- c("Black", "Brown", "Blonde", "Red", "Gray")
  colnames(errorMat) <- c("Black", "Brown", "Blonde", "Red", "Gray")

  return(errorMat)
}
