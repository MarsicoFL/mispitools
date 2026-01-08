#' Get Allele Frequencies in pedtools Format
#'
#' @description
#' Converts allele frequency data from a data frame to a list format
#' compatible with the \pkg{pedtools} package for pedigree analysis.
#'
#' @param region A data frame containing allele frequencies. The first column
#'   should be "Allele" with allele designations, and subsequent columns
#'   should be marker names with frequency values. Available built-in
#'   databases: \code{\link{Argentina}}, \code{\link{Asia}},
#'   \code{\link{Austria}}, \code{\link{BosniaHerz}}, \code{\link{China}},
#'   \code{\link{Europe}}, \code{\link{Japan}}, \code{\link{USA}}.
#'
#' @return A named list where each element represents a genetic marker.
#'   Each marker element is a named numeric vector with allele names and
#'   their corresponding frequencies. This format is directly compatible
#'   with \code{pedtools::setMarkers()}.
#'
#' @details
#' The function transforms the data frame format (rows = alleles, columns =
#' markers) into the list format required by pedtools (one element per marker,
#' named by allele). This enables seamless integration with pedigree
#' likelihood calculations.
#'
#' @seealso
#' \code{\link{Argentina}}, \code{\link{Europe}}, \code{\link{USA}} for
#' available frequency databases, \code{\link{sim_lr_genetic}} for using
#' these frequencies in simulations.
#'
#' @references
#' Marino M, et al. (2009). "Allele frequencies of 15 STRs in an Argentine
#' population sample." Forensic Science International: Genetics Supplement
#' Series. \doi{10.1016/j.fsigss.2009.08.178}
#'
#' @source
#' \doi{10.1016/j.fsigss.2009.08.178}
#' \doi{10.1016/j.fsigen.2016.06.008}
#' \doi{10.1016/j.fsigen.2018.07.013}
#'
#' @export
#' @examples
#' # Convert Argentina database to pedtools format
#' freqs <- get_allele_freqs(Argentina)
#'
#' # Check available markers
#' names(freqs)
#'
#' # Use with pedtools
#' library(pedtools)
#' library(forrel)
#' x <- linearPed(2)
#' x <- setMarkers(x, locusAttributes = freqs[1:5])

get_allele_freqs <- function(region) {

  Freqs <- as.list(region)

  # Set allele names for each marker
  for (i in 2:length(Freqs)) {
    names(Freqs[[i]]) <- Freqs[[1]]
  }

  # Remove the Allele column
  Freqs$Allele <- NULL

  return(Freqs)
}
