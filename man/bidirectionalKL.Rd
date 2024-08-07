% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bidirectionalKL.R
\name{bidirectionalKL}
\alias{bidirectionalKL}
\title{Kullback-Leibler Divergence Calculation for Genetic Markers}
\usage{
bidirectionalKL(data1, data2, minFreq = 1e-10)
}
\arguments{
\item{data1}{DataFrame with allele frequencies for the first population.}

\item{data2}{DataFrame with allele frequencies for the second population.}

\item{minFreq}{Minimum frequency to be considered for unobserved or poorly observed alleles.}
}
\value{
A list containing the Kullback-Leibler divergence from data1 to data2 and vice versa.
}
\description{
This function calculates the Kullback-Leibler divergence for shared genetic markers
between two populations, considering allele frequencies. It normalizes data, adjusts zero frequencies,
and calculates divergence in both directions.
}
\examples{
bidirectionalKL(Argentina, BosniaHerz)
}
