% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EFGLdata.R
\name{readInData}
\alias{readInData}
\title{read in data from a Progeny-style output file
or matrix/dataframe/tibble adn create an EFGLdata object}
\usage{
readInData(
  input,
  genotypeStart = NULL,
  pedigreeColumn = 1,
  nameColumn = 2,
  convertNames = TRUE,
  convertMetaDataNames = TRUE,
  missingAlleles = c("0", "00", "000"),
  guess_max = 10000
)
}
\arguments{
\item{input}{Either a character string to the tab-separated input file with
a header row or a
matrix/dataframe/tibble. Structure of the input: one row per individual.
One column giving
pedigree (population) names, one column giving individual names,
optional additional metadata columns, genotype columns. Pedigree and
individual name columns can be anywhere, if specified. Genotype columns
MUST be consecutive and be the right most columns. Genotypes are given as
two columns per call (diploidy assumed).}

\item{genotypeStart}{The column number that genotypes start at. If not specified,
the first column with a column name ending in ".A1", ".a1", "-A1", or "-a1" is chosen.}

\item{pedigreeColumn}{The column number that contains pedigree (population) names.}

\item{nameColumn}{The column number that contains individual names. These MUST be unique.}

\item{convertNames}{TRUE to convert genotype and pedigree names in the same way
that IDFGEN does (remove special
characters from both and remove "." from genotype names).}

\item{convertMetaDataNames}{TRUE to remove special characters and spaces from metadata
column names. This makes accessing them easier.}

\item{missingAlleles}{a vector of values (not NA) to treat as missing alleles. They will be converted to NA.}

\item{guess_max}{If input is a character string, this is the maximum number of lines to use when
guessing input data types. Making this smaller results in quicker loading, making it larger
can fix some parsing errors}
}
\value{
An EFGLdata object, which is just a list with two elements. The first element
is a tibble with genotype data, the second is a tibble with metadata
}
\description{
read in data from a Progeny-style output file
or matrix/dataframe/tibble adn create an EFGLdata object
}
