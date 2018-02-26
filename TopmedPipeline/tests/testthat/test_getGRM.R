context("getGRM tests")
library(GENESIS)
library(SeqArray)
library(SNPRelate)

.testConfig <- function(type) {
    config <- character()
    if (type == "pcrelate") {
        config["pcrelate_file"] <- system.file("extdata", "HapMap_ASW_MXL_pcrelate.gds", package="GENESIS")
    } else if (type == "grm") {
        gds <- seqOpen(seqExampleFileName())
        grm <- snpgdsGRM(gds, verbose=FALSE)
        seqClose(gds)
        grmfile <- tempfile()
        save(grm, file=grmfile)
        config["grm_file"] <- grmfile
    } else if (type == "grm_gds") {
        gds <- seqOpen(seqExampleFileName())
        grmfile <- file.path(tempdir(), "tmp.gds")
        grm <- snpgdsGRM(gds, out.fn=grmfile, verbose=FALSE)
        seqClose(gds)
        config["grm_file"] <- grmfile
    }
    config
}

.cleanupConfig <- function(config) {
    unlink(config[c("grm_file")])
}

test_that("pcrelate", {
    config <- .testConfig(type="pcrelate")
    pcr <- openfn.gds(config["pcrelate_file"])
    samp <- as.character(read.gdsn(index.gdsn(pcr, "sample.id"))[1:10])
    closefn.gds(pcr)
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "matrix")
    expect_equal(colnames(grm), samp)
})

test_that("grm", {
    config <- .testConfig(type="grm")
    
    x <- getobj(config["grm_file"])
    samp <- x$sample.id[1:10]
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})

test_that("grm_gds", {
    config <- .testConfig(type="grm_gds")
    
    x <- openfn.gds(config["grm_file"])
    samp <- read.gdsn(index.gdsn(x, "sample.id"))[1:10]
    closefn.gds(x)
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})
