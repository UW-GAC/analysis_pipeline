
#' @importFrom dplyr "%>%" arrange_ select_
.arrange_chr_pos <- function(df, chr="chr", pos="pos") {
    df$chr.factor <- factor(df[[chr]], levels=c(1:22, "X", "Y"))
    df %>%
        arrange_("chr.factor", pos) %>%
        select_("-chr.factor")
}

.index_chr_pos <- function(x, chr="chr", pos="pos") {
    df <- do.call(rbind, lapply(x, function(xx) xx[1,c(chr,pos)]))
    df$chr.factor <- factor(df[[chr]], levels=c(1:22, "X", "Y"))
    order(df$chr.factor, df[[pos]])
}


#' Combine association test results
#'
#' Combine association test results from multiple files into a single object.
#'
#' Useful for combining per-segment results into a single file per chromosome.
#' Assumes files are already ordered by segment.
#'
#' @param files Vector of file names with association test results
#' @param assoc_type Type of association test ("single", "aggregate", "window")
#' @param ordered Logical for whether to order the output by chromosome and position
#' @return Association test object
#'
#' @importFrom dplyr "%>%" bind_rows distinct_ filter_ group_by_ mutate_
#' @export
combineAssoc <- function(files, assoc_type, ordered=FALSE) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    x <- lapply(unname(files), getobj)
    if (assoc_type == "single") {
        assoc <- bind_rows(x)
        if (ordered) {
            assoc <- .arrange_chr_pos(assoc)
        }
    } else if (assoc_type %in% c("aggregate", "window")) {
        assoc <- list()
        #assoc <- x[[1]][c("param", "nsample")]
        assoc$results <- bind_rows(lapply(x, function(y) y$results))
        assoc$variantInfo <- do.call(c, lapply(x, function(y) y$variantInfo))
        if (ordered) {
            # get index to put units in order by chr, pos
            index <- .index_chr_pos(assoc$variantInfo)
            assoc$results <- assoc$results[index,]
            assoc$variantInfo <- assoc$variantInfo[index]
        }
    }
    if (assoc_type == "aggregate") {
        if (!is.null(names(assoc$variantInfo))) {
            # remove duplicated units
            index <- !duplicated(names(assoc$variantInfo))
            assoc$results <- assoc$results[index,]
            assoc$variantInfo <- assoc$variantInfo[index]
        }
    }
    if (assoc_type == "window") {
        # remove duplicated windows
        assoc$results <- assoc$results %>%
            mutate_(index=~1:n()) %>%
            group_by_("chr", "start", "end") %>%
            filter_(~(n.site == max(n.site)), ~(!duplicated(n.site))) %>%
            as.data.frame()
        assoc$variantInfo <- assoc$variantInfo[assoc$results$index]
        assoc$results <- assoc$results %>%
            select_(~-index)
    }
    assoc
}

#' Get association test results
#'
#' Return association test results in a standard format
#'
#' Read association test results in multiple files and combine all into a single
#' data frame with standard column names.
#'
#' The \code{id} values in this file should be either:
#' \itemize{
#'   \item{single-variant tests: }{\code{variant.id}}
#'   \item{aggregate tests: }{the name of the aggregation unit from the variant grouping file}
#'   \item{window tests: }{\code{<chr>_<pos>} of the first variant, e.g., 2_20001 for a window on chr2 starting at position 200001}
#' }
#'
#' If a single aggregate unit contains variants from multiple chromosomes, the result will be returned only once on the chromosome with the most variants.
#'
#' @inheritParams combineAssoc
#' @return data.frame including standard columns ("id", "chr", "pos", "start", "end", "stat", "pval", "MAC"). Also includes "MAF" for single variant tests.
#'
#' @importFrom dplyr "%>%" bind_rows ends_with filter group_by left_join mutate n rename select summarise .data arrange desc
#' @export
getAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    assoc <- bind_rows(lapply(unname(files), function(f) {
        x <- getobj(f)
        if (assoc_type == "aggregate") {
            tmp <- x$results %>%
                mutate(group_id=names(x$variantInfo)) %>%
                filter(.data$n.site > 0)
            group.info <- x$variantInfo %>%
              bind_rows(.id="group_id") %>%
              group_by(.data$group_id, .data$chr) %>%
              summarise(n_chr = n(),
                        start=min(.data$pos),
                        end=max(.data$pos),
                        pos=(floor((min(.data$pos) + max(.data$pos))/2))) %>%
              ungroup() %>%
              group_by(.data$group_id) %>%
              # Choose the chromosome with the most variants as the one to return.
              arrange(desc(.data$n_chr)) %>%
              dplyr::slice(1) %>%
              ungroup()
            x <- left_join(tmp, group.info, by="group_id") %>%
              mutate(id=.data$group_id)
        } else if (assoc_type == "window") {
            x <- filter(x$results, (.data$n.site > 0)) %>%
                mutate(pos=(floor((.data$start + .data$end)/2))) %>%
                mutate(id = sprintf("chr%s_%d", .data$chr, .data$start))
        } else {
            x <- mutate(x, start=.data$pos, end=.data$pos) %>%
              dplyr::rename(id = .data$variant.id)
        }
        x
    }))

    # This could probably be cleaned up to look at assoc_type first and then select specific p-value.
    if (assoc_type == "single" && "mid.pval" %in% names(assoc)) {
      ## BinomiRare
      # Select the mid.pval as recommended by Tamar.
      # No test stat.
      assoc <- assoc %>%
        mutate(MAF = pmin(1 - freq, freq)) %>%
        select(.data$id, .data$chr, .data$pos, .data$start, .data$end, .data$mid.pval, .data$MAC, .data$MAF)
      names(assoc)[6] <- c("pval")
    } else if ("pval" %in% names(assoc)) {
        ## SKAT or fastSKAT
        assoc <- select(assoc, .data$id, .data$chr, .data$pos, .data$start, .data$end, .data$pval, suppressWarnings(one_of("MAC")))
    } else if ("pval_SKATO" %in% names(assoc)) {
        ## SKATO
        assoc <- select(assoc, .data$id, .data$chr, .data$pos, .data$start, .data$end, pval=.data$pval_SKATO, suppressWarnings(one_of("MAC")))
    } else if ("pval_SMMAT" %in% names(assoc)) {
        ## SMMAT
        assoc <- select(assoc, .data$id, .data$chr, .data$pos, .data$start, .data$end, pval=.data$pval_SMMAT, suppressWarnings(one_of("MAC")))
    } else if (assoc_type == "single") {
      assoc <- assoc %>%
        mutate(MAF = pmin(1 - freq, freq)) %>%
        select(.data$id, .data$chr, .data$pos, .data$start, .data$end, ends_with("Stat"), ends_with("pval"), .data$MAC, .data$MAF)
      names(assoc)[6:7] <- c("stat", "pval")
    } else {
        ## burden
        assoc <- select(assoc, .data$id, .data$chr, .data$pos, .data$start, .data$end, ends_with("Stat"), ends_with("pval"),
                         suppressWarnings(one_of("MAC")))
        names(assoc)[6:7] <- c("stat", "pval")
    }

    assoc <- filter(assoc, !is.na(.data$pval)) %>%
        mutate(chr=ordered(.data$chr, levels=c(1:22, "X")))
    assoc
}


#' Omit known hits from an association test data frame
#'
#' @param assoc data.frame with assocation test results (including columns chr, pos)
#' @param hits data.frame with known hits (including columns chr, pos)
#' @param flank Number of kb on either side of each known hit to exclude
#' @return data.frame with assocation test results not in regions around known hits
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @export
omitKnownHits <- function(assoc, hits, flank=500) {
    stopifnot(all(c("chr", "pos") %in% names(hits)))
    assoc.gr <- GRanges(seqnames=assoc$chr, ranges=IRanges(start=assoc$pos, end=assoc$pos))
    hits.gr <- GRanges(seqnames=hits$chr, ranges=IRanges(start=hits$pos-(flank*1000),
                                              end=hits$pos+(flank*1000)))
    ol <- findOverlaps(assoc.gr, hits.gr)
    assoc[-queryHits(ol),]
}


#' Add MAC column to association test output
#'
#' @param assoc results from \code{\link[GENESIS]{assocTestSingle}} or \code{\link[GENESIS]{assocTestAggregate}}
#' @param assoc_type Type of association test ("single", "aggregate", "window")
#' @return \code{assoc} with a "MAC" column added to the results data.frame
#'
#' @export
addMAC <- function(assoc, assoc_type) {
    mac <- function(x) {
        if ("MAC" %in% names(x)) {
            x$MAC
        } else {
            round(2 * x$n.obs * pmin(x$freq, 1-x$freq))
        }
    }
    if (assoc_type == "single") {
        assoc$MAC <- mac(assoc)
    } else if (assoc_type %in% c("aggregate", "window")) {
        assoc$results$MAC <- sapply(assoc$variantInfo, function(x) sum(mac(x)))
    }
    assoc
}


#' Remove conditional variants from assoc file
#'
#' @param assoc results from \code{\link[GENESIS]{assocTestSingle}}
#' @param varfile conditional variant filename (RData containing data.frame with columns "variant.id", "chr")
#'
#' @return \code{assoc} with conditional variants removed
#' @importFrom dplyr anti_join select_
#'
#' @export
removeConditional <- function(assoc, varfile) {
    dat <- getobj(varfile)
    if ("chromosome" %in% names(dat)) names(dat)[names(dat) == "chromosome"] <- "chr"
    stopifnot(all(c("chr", "variant.id") %in% names(dat)))
    dat <- select_(dat, "variant.id", "chr")
    dat$chr <- as.character(dat$chr)
    anti_join(assoc, dat, by=c("variant.id", "chr"))
}

#' Filter association test results to just the ids specified in a file.
#' @param assoc results from \code{\link{getAssoc}}
#' @param varfile RData file containing vector of \code{id} values from \code{assoc} to keep. See details.
#'
#' @return \code{assoc} with only the specified identifiers
#'
#' @details
#' If \code{varfile} contains a space, `chromosome` will be inserted and the variant ids included will be assumed to be from that chromosome.
#' Ids should match the \code{id} column returned by \code{\link{getAssoc}}.
#'
#' @importFrom dplyr inner_join mutate filter
#'
#' @export
#'
assocFilterByFile <- function(assoc, varfile) {
  # Is the variant include file by chromosome?
  if (grepl(" ", varfile)) {
    chrs <- as.character(unique(assoc$chr))
    keep <- lapply(chrs, function(x) {
      varfile_chr <- insertChromString(varfile, x)
      if (file.exists(varfile_chr)) {
        tmp <- getobj(varfile_chr)
      } else {
        warning(sprintf("missing varfile for chromosome %s; no variants included from this chromosome", x))
        tmp <- NULL
      }
      # Create an empty vector for cases when there are no variants specified in this file.
      if (is.null(tmp)) tmp <- get(class(assoc$id))()
      data.frame(id = tmp, chr = rep(x, length(tmp)), stringsAsFactors = FALSE)
    }) %>%
      bind_rows()
    keep <- keep %>%
      # Convert to factor ordering for the chromosome.
      mutate(chr = ordered(chr, levels = levels(assoc$chr)))
    assoc <- assoc %>%
      inner_join(keep, by = c("id", "chr"))
  } else {
    var_include <- getobj(varfile)
    assoc <- assoc %>%
      filter(id %in% var_include)
  }

  assoc
}
