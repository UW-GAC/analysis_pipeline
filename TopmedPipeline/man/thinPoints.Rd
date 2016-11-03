\name{thinPoints}
\alias{thinPoints}
\title{Thin points in a data frame for plotting}
\description{Sample a data frame by sampling the value of one of the columns.}
\usage{thinPoints(dat, value, n=10000, nbins=10, groupBy=NULL)}
\arguments{
  \item{dat}{a data.frame}
  \item{value}{the column in \code{dat} to sample on}
  \item{n}{the number of points to select in each bin}
  \item{nbins}{the number of bins to divide \code{dat$value} into}
  \item{groupBy}{if not \code{NULL}, bins are calculated within each value of \code{groupBy} separately}
}
%\details{}
\value{A subset of \code{dat}. The number of rows will equal \code{n*nbins*ngroup} the number of unique elements in \code{dat$groupBy}.}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
dat <- data.frame(x=1:100, y=sample(letters[1:2], 100, replace=TRUE))
thinPoints(dat, "x", n=2, nbins=5, groupBy="y")
}
\keyword{manip}
