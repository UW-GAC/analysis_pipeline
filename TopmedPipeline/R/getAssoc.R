getAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    do.call(rbind, lapply(unname(files), function(f) {
        assoc <- getobj(f)
        if (assoc_type  == "aggregate") {
            tmp <- assoc$results %>%
                rownames_to_column("group_id") %>%
                filter_(~(n.site > 0))
            group.info <- do.call(rbind, lapply(tmp$group_id, function(g) {
                grp <- assoc$variantInfo[[g]][1, c("chr", "pos")]
                grp$group_id <- g
                grp
            }))
            assoc <- left_join(tmp, group.info, by="group_id")
        } else if (assoc_type == "window") {
            assoc <- filter_(assoc$results, ~(n.site > 0), ~(dup == 0)) %>%
                mutate_(pos=~(floor((window.start + window.stop)/2)))
        }
        assoc
    }))
}
