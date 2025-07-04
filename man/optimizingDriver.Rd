% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Driver.R
\name{optimizingDriver}
\alias{optimizingDriver}
\title{Run \code{reverseGWAS} on an input dataset in full optimisation mode for a single genotype, i.e. to find the best association to within a given factor.}
\usage{
optimizingDriver(
  inputFile,
  startPValue = log(0.001),
  factor = 2,
  type = "CNF",
  complement = -1,
  K = MY_K,
  L = MY_L,
  startIndex = 0
)
}
\arguments{
\item{inputFile}{The absolute path to an input file. This can be a CSV or TSV file, either uncompressed or compressed with gz. It must start with a patient identifier column named ID (typically a character), followed by one or more genotype columns, followed by one or more phenotype columns. NB: the file extension must be immediately preceded by an integer stating the number of genotype columns.}

\item{startPValue}{Numeric value, the logarithm of a minimum acceptable p-value (if infeasible, no further optimisation is carried out); the default is \code{log(1e-3)}.}

\item{factor}{Numeric value, a multiplicative factor to within which the best association is sought; the default is \code{2}.}

\item{type}{The type of phenotype combinations to search over; the allowed values are "CNF", or conjunctive normal form (the default) and "DNF", or disjunctive normal form.}

\item{complement}{Integer between \code{-2} and \code{1}, or a value coercible to one. If it is  \code{1} or \code{2}, each genotype's negated form is added in preprocessing. If it is \code{-1} or \code{2}, each phenotype's negated form is added in preprocessing.}

\item{K, L}{The maximum number of clauses (K) and the maximum size of each clause (L) to use in searching for the combination phenotype. The defaults are \code{2} and \code{2}. NB: high values may lead to overfitting to the training (e.g. discovery) data!}

\item{startIndex}{Offset value to be used in distinguishing between the iterations; the default is \code{0}.}
}
\value{
A list containing a summary of the association statistics (one row per genotype) and the best combination phenotypes (one column per genotype).
}
\description{
Carry out the discovery-validation pipeline, whose arguments parallel those of the mainDriver and similar functions.
}
\examples{
# Z <- optimizingDriver(inputFile = "TestInputN5000P5_1.csv", startPValue = log(1e-3))
}
