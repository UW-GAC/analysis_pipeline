library(argparser)
library(GWASTools)
library(dplyr)
library(tidyr)
sessionInfo()

argp <- arg_parser("Pedigree check")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "pedigree_file", help="path to pedigree file")
argv <- parse_args(argp)
config <- readConfig(argv$config)

required <- c("pedigree_file")
optional <- c("out_file"="exp_rel.RData",
              "err_file"="ped_errs.RData")
config <- setConfigDefaults(config, required, optional)
print(config)

config <- c(pedigree_file=argv$pedigree_file,
            out_file="out.RData", err_file="err.RData")

pedfile <- config["pedigree_file"]

sep <- switch(tools::file_ext(pedfile),
              csv = ",",
              txt = "\t",
              "")
hdr <- !(tools::file_ext(pedfile) == "fam")
ped <- read.table(pedfile, sep=sep, header=hdr, comment.char="#", na.strings=c("", "NA"), as.is=TRUE, fill=TRUE, blank.lines.skip=TRUE)
head(ped)

# if this is a dbGaP file, strip out dbGaP subject ID so we only have one ID column
if (any(grepl("^dbGaP", names(ped)))) {
    ped <- ped %>%
        select(-starts_with("dbGaP"))
}

# check column names
names(ped) <- tolower(names(ped))
cols <- lapply(c("fam", "subj", "father", "mother", "sex"),
               function(x) {
                   i <- which(grepl(x, names(ped)))
                   if (length(i) == 1) return(i) else return(NA)
               })
names(cols) <- c("family", "individ", "father", "mother", "sex")
cols <- unlist(cols)

if (is.na(cols["family"])) cols["family"] <- 1
if (is.na(cols["individ"])) cols["individ"] <- 2
if (is.na(cols["father"])) cols["father"] <- 3
if (is.na(cols["mother"])) cols["mother"] <- 4
if (is.na(cols["sex"])) cols["sex"] <- 5

if (!setequal(cols, 1:5)) {
    stop("Cannot parse pedigree file. Columns should be FAMILY_ID, SUBJECT_ID, FATHER, MOTHER, SEX.")
}

ped <- ped[,unname(cols)]
names(ped) <- names(cols)


# set mother and father ID for founders to 0
ped <- ped %>%
    mutate(father=ifelse(is.na(father), 0, father),
           mother=ifelse(is.na(mother), 0, mother))


# standardize sex column
if (is.numeric(ped$sex)) {
    ped <- ped %>%
        mutate(sex=c("M", "F")[sex])
} else {
    if (any(tolower(ped$sex) %in% c("male", "female"))) {
        ped <- ped %>%
            mutate(sex=toupper(substr(sex,1,1)))
    }
}
head(ped)

# make sure individ is unique
if (any(duplicated(ped$individ))) {
    ped <- ped %>%
        mutate(individ=paste(family, individ, sep="_"),
               father=ifelse(father == 0, 0, paste(family, father, sep="_")),
               mother=ifelse(mother == 0, 0, paste(family, mother, sep="_")))
}


# check for pedigree errors
chk <- pedigreeCheck(ped)
names(chk)

if (!is.null(chk)) {
    save(chk, file=config["err_file"])
}

if ("duplicates" %in% names(chk)) {
    if (!all(chk$duplicates$match)) {
        stop("Pedigree has duplicate subjects with conflicting data.")
    }
    ped <- distinct(ped)
}

if ("unknown.parent.rows" %in% names(chk)) {
    ped <- ped %>%
        mutate(father=ifelse(father == 0 & mother != 0,
                             paste(family, "father", row_number(), sep="_"),
                             father),
               mother=ifelse(father != 0 & mother == 0,
                             paste(family, "mother", row_number(), sep="_"),
                             mother)) 
}


## repeat check after correcting unknown parent rows to get new
## dummy parents in "no.individ.entry"
chk <- pedigreeCheck(ped)

if ("parent.no.individ.entry" %in% names(chk)) {
    both <- chk$parent.no.individ.entry %>%
        filter(no_individ_entry == "both") %>%
        separate(parentID, into=c("mother", "father"), sep=";") %>%
        select(family, father, mother) %>%
        pivot_longer(-family, names_to="no_individ_entry", values_to="parentID")
    parents <- chk$parent.no.individ.entry %>%
        filter(no_individ_entry != "both") %>%
        bind_rows(both) %>%
        mutate(sex=ifelse(no_individ_entry == "father", "M", "F"),
               father="0", mother="0") %>%
        select(family, individ=parentID, father, mother, sex)
    ped <- bind_rows(ped, parents)
}

## repeat check after adding dummy parents
chk <- pedigreeCheck(ped)

if ("one.person.fams" %in% names(chk)) {
    ped <- ped %>%
        filter(!family %in% chk$one.person.fams$family)
}

if ("subfamilies.ident" %in% names(chk)) {
    ped <- ped %>%
        left_join(chk$subfamilies.ident, by=c("family", "individ")) %>%
        mutate(family=ifelse(is.na(subfamily), family, paste(family, subfamily, sep="_"))) %>%
        select(-subfamily)
}

chk <- pedigreeCheck(ped)

## sometimes we need to do this again after assigning subfamilies
if ("one.person.fams" %in% names(chk)) {
    ped <- ped %>%
        filter(!family %in% chk$one.person.fams$family)
}

chk <- pedigreeCheck(ped)
names(chk)


if (!is.null(chk)) {
    stop("pedigree had unresolvable errors")
}


# define relative categories
source("https://raw.githubusercontent.com/UW-GAC/QCpipeline/master/QCpipeline/R/expRelsCategory.R")

rel <- expRelsCategory(ped)
rel <- rel$relprs.all

save(rel, file=config["out_file"])

table(rel$relation, rel$exp.rel)
