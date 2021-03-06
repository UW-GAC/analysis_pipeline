% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAssoc.R
\name{getAssoc}
\alias{getAssoc}
\title{Get association test results}
\usage{
getAssoc(files, assoc_type)
}
\arguments{
\item{files}{Vector of file names with association test results}

\item{assoc_type}{Type of association test ("single", "aggregate", "window")}
}
\value{
data.frame including standard columns ("id", "chr", "pos", "start", "end", "stat", "pval", "MAC"). Also includes "MAF" for single variant tests.
For BinomiRare tests, the `pval` column in the returned data frame is the `mid.pval`.
}
\description{
Return association test results in a standard format
}
\details{
Read association test results in multiple files and combine all into a single
data frame with standard column names.

The \code{id} values in this file should be either:
\itemize{
  \item{single-variant tests: }{\code{variant.id}}
  \item{aggregate tests: }{the name of the aggregation unit from the variant grouping file}
  \item{window tests: }{\code{<chr>_<pos>} of the first variant, e.g., 2_20001 for a window on chr2 starting at position 200001}
}

If a single aggregate unit contains variants from multiple chromosomes, the result will be returned only once on the chromosome with the most variants.
}
