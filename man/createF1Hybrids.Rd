% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{createF1Hybrids}
\alias{createF1Hybrids}
\title{Create F1 hybrids from two populations}
\usage{
createF1Hybrids(x, pop1, pop2, newName, n = 50, missingGenos = TRUE)
}
\arguments{
\item{x}{an EFGLdata object}

\item{pop1}{one parent population (must be in `x`)}

\item{pop2}{the other parent population (must be in `x`)}

\item{newName}{a string giving the name of the population to put hybrids
into. This MUST be a new pop.}

\item{n}{the number of hybrids to simulate}

\item{missingGenos}{TRUE to simulate missing genotypes, FALSE to not}
}
\value{
an EFGLdata object with all the populations in `x`, plus the new simulated hybrid population
}
\description{
Creates hybrids by sampling one allele from one pop
and the other allele from the other pop
based on observed allele frequencies. Missing genotypes
are sampled based on the combined rate of missing genotypes
or data can be simulated wiht no missing genotypes. Alleles
at different loci are sampled independently. Note that if a sex marker exists
in the data set, it is treated like any other marker, resulting in XX, XY, and YY
genotypes in the hybrids. You may want to remove it prior to simulation using
the `removeLoci` function. Loci with all missing genotypes in one population
are simulated to have all missing genotypes in the hybrids REGARDLESS of
the input value for `missingGenos`.
}
