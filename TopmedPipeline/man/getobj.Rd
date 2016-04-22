\name{getobj}
\alias{getobj}
\title{Get an R object stored in an Rdata file}
\description{Returns an R object stored in an Rdata file}
\usage{getobj(Rdata)}
\arguments{
  \item{Rdata}{path to an Rdata file containing a single R object to load}
}
\details{Loads an R object and stores it under a new name without
  creating a duplicate copy.  If multiple objects are stored in the same
  file, only the first one will be returned}
\value{The R object stored in \code{Rdata}.}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
x <- 1:10
file <- tempfile()
save(x, file=file)
y <- getobj(file)
unlink(file)
}
\keyword{utilities}
