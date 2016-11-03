## thin points in a data frame for plotting
## http://stackoverflow.com/questions/30950016/dplyr-sample-n-where-n-is-the-value-of-a-grouped-variable

thinPoints <- function(dat, value, n=10000, nbins=10, groupBy=NULL){
    if (!is.null(groupBy)) {
        dat <- group_by_(dat, groupBy)
    }

    dat %>%
        mutate_(bin=interp(~cut(value, breaks=nbins, labels=FALSE), value=as.name(value), nbins=nbins)) %>%
        group_by_(~bin, add=TRUE) %>%
        sample_frac(1) %>%
        filter_(~(row_number() <= min(n)), ~n()) %>%
        ungroup() %>%
        select_(~(-bin))
}
