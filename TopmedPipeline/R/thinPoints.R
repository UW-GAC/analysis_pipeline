#' Thin points in a data frame for plotting
#'
#' Sample a data frame by sampling the value of one of the columns
#'
#' @param dat a data.frame
#' @param value the column in \code{dat} to sample on
#' @param n the number of points to select in each bin
#' @param nbins the number of bins to divide \code{dat$value} into
#' @param groupBy if not \code{NULL}, bins are calculated within each value of \code{groupBy} separately
#' @return A subset of \code{dat}. The maxumum number of rows is \code{n*nbins*ngroup} where
#'   \code{ngroup} = the number of unique elements in \code{dat$groupBy}. The number of rows may be less
#'   if some bins have fewer than \code{n} elements.
#' @examples
#' dat <- data.frame(x=1:100, y=sample(letters[1:2], 100, replace=TRUE))
#' thinPoints(dat, "x", n=2, nbins=5, groupBy="y")
#' @references \url{http://stackoverflow.com/questions/30950016/dplyr-sample-n-where-n-is-the-value-of-a-grouped-variable}
#'
#' @importFrom dplyr "%>%" filter_ group_by_ mutate_ row_number sample_frac select_ ungroup
#' @importFrom lazyeval interp
#' @export
thinPoints <- function(dat, value, n=10000, nbins=10, groupBy=NULL){
    if (!is.null(groupBy)) {
        dat <- group_by_(dat, groupBy)
    }

    dat %>%
        mutate_(bin=interp(~cut(value, breaks=nbins, labels=FALSE), value=as.name(value), nbins=nbins)) %>%
        group_by_(~bin, add=TRUE) %>%
        sample_frac(1) %>%
        filter_(~(row_number() <= n)) %>%
        ungroup() %>%
        select_(.dots="-bin")
}
