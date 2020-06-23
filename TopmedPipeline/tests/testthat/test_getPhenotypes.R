context("getPhenotype tests")
library(splines)
library(Biobase)
library(SeqArray)
library(SeqVarTools)
library(gdsfmt)

.testConfig <- function(covars, n_pcs, group_var=NA) {
    n <- 100
    phen <- AnnotatedDataFrame(data.frame(sample.id=as.character(1:n),
                                          outcome=rnorm(n),
                                          study=sample(letters, n, replace=TRUE),
                                          age=runif(n,18,70),
                                          stringsAsFactors=FALSE))
    phen$age.spline <- bs(phen$age, df=3)
    phenfile <- tempfile()
    save(phen, file=phenfile)

    pcs <- list(vectors=matrix(rnorm(n*2), nrow=n, ncol=2, dimnames=list(phen$sample.id), 1:2))
    class(pcs) <- "pcair"
    pcfile <- tempfile()
    save(pcs, file=pcfile)

    c(phenotype_file=phenfile,
      pca_file=pcfile,
      n_pcs=n_pcs,
      outcome="outcome",
      covars=covars,
      group_var=group_var,
      sample_include_file=NA)
}

.cleanupConfig <- function(config) {
    unlink(config[c("phenotype_file", "pca_file")])
}

test_that("one covariate", {
    config <- .testConfig(covars="study", n_pcs=2)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_equal(phen$covars, c("study", "PC1" ,"PC2"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study", "PC1" ,"PC2"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("no PCs", {
    config <- .testConfig(covars="study", n_pcs=0)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_equal(phen$covars, c("study"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("no covariates", {
    config <- .testConfig(covars=NA, n_pcs=2)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_equal(phen$covars, c("PC1", "PC2"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "PC1", "PC2"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("no covariates or PCs", {
    config <- .testConfig(covars=NA, n_pcs=0)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_null(phen$covars)
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("no PC file error", {
    config <- .testConfig(covars=NA, n_pcs=2)
    config["pca_file"] <- NA
    expect_error(getPhenotypes(config), "must specify pca_file")
    
    .cleanupConfig(config)
})

test_that("one PC", {
    config <- .testConfig(covars=NA, n_pcs=1)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_equal(phen$covars, "PC1")
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "PC1"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("group.var is NA", {
    config <- .testConfig(covars="study", n_pcs=0, group_var=NA)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$covars, c("study"))
    expect_null(phen$group.var)
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("group.var is covar", {
    config <- .testConfig(covars="study", n_pcs=0, group_var="study")
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$covars, c("study"))
    expect_equal(phen$group.var, c("study"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("group.var is not covar", {
    config <- .testConfig(covars=NA, n_pcs=0, group_var="study")
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_null(phen$covars)
    expect_equal(phen$group.var, c("study"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})


test_that("spline", {
    config <- .testConfig(covars="study age.spline", n_pcs=2)
    
    phen <- getPhenotypes(config)
    expect_is(phen$annot, "AnnotatedDataFrame")
    expect_equal(phen$outcome, "outcome")
    expect_equal(phen$covars, c("study", "age.spline", "PC1" ,"PC2"))
    expect_equal(varLabels(phen$annot), c("sample.id", "outcome", "study", "age.spline", "PC1", "PC2"))
    expect_equal(phen$sample.id, phen$annot$sample.id)

    .cleanupConfig(config)
})

test_that("fewer samples with pcs", {
    config <- .testConfig(covars="study", n_pcs=2)
    pcs <- getobj(config["pca_file"])
    keep <- rownames(pcs$vectors)[sort(sample.int(100,50))]
    pcs$vectors <- pcs$vectors[keep,]
    save(pcs, file=config["pca_file"])
    
    phen <- getPhenotypes(config)
    expect_equal(phen$sample.id, keep)
    expect_true(all(is.na(phen$annot$PC1[!(phen$annot$sample.id %in% keep)])))

    .cleanupConfig(config)
})

test_that("parseParam", {
    expect_null(.parseParam(NA))
    expect_equal(.parseParam("a"), "a")
    expect_equal(.parseParam("a b c"), c("a", "b", "c"))
})

test_that("genotypes", {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    variant.id <- sample(1:100, 2)
    geno <- .genotypes(gdsfile, variant.id)
    expect_true(is(geno, "matrix"))
    expect_equal(colnames(geno), as.character(sort(variant.id)))
})

.testConditionalConfig <- function(gdsfile, nchr=1, nvar=3) {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    gds <- seqOpen(gdsfile)
    dat <- data.frame(variant.id=seqGetData(gds, "variant.id"),
                      chromosome=seqGetData(gds, "chromosome"),
                      stringsAsFactors=FALSE)
    seqClose(gds)
    chr <- sample(unique(dat$chromosome), nchr)
    dat <- dat[dat$chromosome %in% chr,]
    dat <- dat[sample(1:nrow(dat), nvar),]
    cvfile <- tempfile()
    save(dat, file=cvfile)
    if (nchr > 1) {
        gdsfile <-  file.path(tempdir(), "tmp .gds")
        for (c in chr) {
            invisible(file.copy(seqExampleFileName("gds"), paste0(tempdir(), "/tmp", c, ".gds")))
        }
    }
    c(gds_file=gdsfile,
      phenotype_file=system.file("data", "sample_annotation.RData", package="TopmedPipeline"),
      conditional_variant_file=cvfile,
      n_pcs=0,
      outcome="outcome",
      covars="sex",
      sample_include_file=NA)
}

.cleanupConditionalConfig <- function(config) {
    unlink(config["conditional_variant_file"])
    unlink(file.path(tempdir(), "/tmp*.gds"))
}
    

test_that("conditional variants", {
    config <- .testConditionalConfig(seqExampleFileName("gds"), nchr=1, nvar=3)
    phen <- getPhenotypes(config)
    cv <- getobj(config["conditional_variant_file"])
    vars <- paste0("var_", cv$variant.id)
    expect_true(all(vars %in% varLabels(phen$annot)))
    expect_true(all(vars %in% phen$covars))
    .cleanupConditionalConfig(config)
})

test_that("conditional chroms", {
    config <- .testConditionalConfig(seqExampleFileName("gds"), nchr=2, nvar=5)
    phen <- getPhenotypes(config)
    cv <- getobj(config["conditional_variant_file"])
    vars <- paste0("var_", cv$variant.id)
    expect_true(all(vars %in% varLabels(phen$annot)))
    expect_true(all(vars %in% phen$covars))
    .cleanupConditionalConfig(config)
})


.testSmallConfig <- function(gdsfile) {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <-  file.path(tempdir(), "tmp .gds")
    invisible(file.copy(seqExampleFileName("gds"), paste0(tempdir(), "/tmp1.gds")))
    
    phenfile <- tempfile()
    phen <- getobj(system.file("data", "sample_annotation.RData", package="TopmedPipeline"))
    phen <- phen[sort(sample(1:nrow(phen), 50)),]
    save(phen, file=phenfile)

    c(gds_file=gdsfile,
      phenotype_file=phenfile,
      n_pcs=0,
      outcome="outcome",
      covars="sex",
      sample_include_file=NA)
}

.cleanupSmallConfig <- function(config) {
    unlink(config["phenotype_file"])
    unlink(file.path(tempdir(), "/tmp*.gds"))
}
    

test_that("fewer samples in phenotype file", {
    config <- .testSmallConfig(seqExampleFileName("gds"))
    phen <- getPhenotypes(config)
    expect_equal(sum(is.na(phen$annot$outcome)), 40)
    svd <- SeqVarData(insertChromString(config["gds_file"], 1), phen$annot)
    seqClose(svd)
    .cleanupSmallConfig(config)
})


test_that("pca_file is matrix or data.frame", {
    config <- .testConfig(covars=NA, n_pcs=2)
    pcafile <- config["pca_file"]
    pcs <- getobj(pcafile)
    pcs <- pcs$vectors
    save(pcs, file=pcafile)
    phen <- getPhenotypes(config)
    expect_equivalent(as.matrix(pData(phen$annot)[,c("PC1", "PC2")]), pcs)

    pcs <- as.data.frame(pcs)
    save(pcs, file=pcafile)
    phen <- getPhenotypes(config)
    expect_equivalent(pData(phen$annot)[,c("PC1", "PC2")], pcs)
    
    .cleanupConfig(config)
})
