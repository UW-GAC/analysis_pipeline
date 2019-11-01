context("getGRM tests")
library(GENESIS)
library(SeqArray)
library(SNPRelate)
library(Matrix)

.testConfig <- function(type) {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    config <- character()
    if (type == "pcrelate") {
        gds <- seqOpen(seqExampleFileName())
        iter <- SeqVarTools::SeqVarBlockIterator(SeqVarTools::SeqVarData(gds), verbose=FALSE)
        pca <- pcair(iter, verbose=FALSE)
        SeqVarTools::resetIterator(iter, verbose=FALSE)
        grm <- pcrelate(iter, pcs=pca$vectors[,1:2], verbose=FALSE)
        seqClose(gds)
        grm <- pcrelateToMatrix(grm, scaleKin=2, verbose=FALSE)
        grmfile <- tempfile()
        save(grm, file=grmfile)
    } else if (type == "grm") {
        gds <- seqOpen(seqExampleFileName())
        grm <- snpgdsGRM(gds, verbose=FALSE)
        seqClose(gds)
        grmfile <- tempfile()
        save(grm, file=grmfile)
    } else if (type == "grm_gds") {
        gds <- seqOpen(seqExampleFileName())
        grmfile <- file.path(tempdir(), "tmp.gds")
        grm <- snpgdsGRM(gds, out.fn=grmfile, verbose=FALSE)
        seqClose(gds)
    } else if (type == "grm_Matrix") {
        gds <- seqOpen(seqExampleFileName())
        grm <- snpgdsGRM(gds, verbose=FALSE)
        seqClose(gds)
        grm <- Matrix(grm$grm, dimnames=list(grm$sample.id, grm$sample.id))
        grmfile <- tempfile()
        save(grm, file=grmfile)
    }
    config["relatedness_matrix_file"] <- grmfile
    config
}

.cleanupConfig <- function(config) {
    unlink(config[c("relatedness_matrix_file")])
}

test_that("pcrelate", {
    config <- .testConfig(type="pcrelate")
    pcr <- getobj(config["relatedness_matrix_file"])
    samp <- rownames(pcr)[1:10]
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "list")
    grm <- grm[[1]]
    expect_is(grm, "Matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})


test_that("grm", {
    config <- .testConfig(type="grm")
    
    x <- getobj(config["relatedness_matrix_file"])
    samp <- x$sample.id[1:10]
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "list")
    grm <- grm[[1]]
    expect_is(grm, "matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})

test_that("grm_gds", {
    config <- .testConfig(type="grm_gds")
    
    x <- openfn.gds(config["relatedness_matrix_file"])
    samp <- read.gdsn(index.gdsn(x, "sample.id"))[1:10]
    closefn.gds(x)
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "list")
    grm <- grm[[1]]
    expect_is(grm, "matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})

test_that("grm_Matrix", {
    config <- .testConfig(type="grm_Matrix")
    
    x <- getobj(config["relatedness_matrix_file"])
    samp <- rownames(x)[1:10]
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "list")
    grm <- grm[[1]]
    expect_is(grm, "Matrix")
    expect_equal(colnames(grm), samp)
    
    .cleanupConfig(config)
})

test_that("multiple files", {
    config <- .testConfig(type="grm")
    
    x <- getobj(config["relatedness_matrix_file"])
    samp <- x$sample.id[1:10]

    config["relatedness_matrix_file"] <- paste(config["relatedness_matrix_file"], config["relatedness_matrix_file"])
    grm <- getGRM(config, sample.id=samp)
    expect_is(grm, "list")
    expect_equal(length(grm), 2)
    expect_equal(grm[[1]], grm[[2]])
    
    .cleanupConfig(config)
})
