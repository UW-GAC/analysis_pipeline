context("getPhenotype tests")
library(splines)
library(Biobase)
library(SeqArray)
library(gdsfmt)

.testConfig <- function(covars, n_pcs) {
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
    pcfile <- tempfile()
    save(pcs, file=pcfile)

    c(phenotype_file=phenfile,
      pca_file=pcfile,
      n_pcs=n_pcs,
      outcome="outcome",
      covars=covars,
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
    expect_true(is(geno, "data.frame"))
    expect_equal(names(geno), c("sample.id", paste0("var_", sort(variant.id))))
})

test_that("conditional variants", {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    config <- c(gds_file=seqExampleFileName("gds"),
                phenotype_file=system.file("data", "sample_annotation.RData", package="TopmedPipeline"),
                conditional_variants="1 2",
                n_pcs=0,
                outcome="outcome",
                covars="sex",
                sample_include_file=NA)
    phen <- getPhenotypes(config)
    expect_true(all(c("var_1", "var_2") %in% varLabels(phen$annot)))
    expect_true(all(c("var_1", "var_2") %in% phen$covars))
})

test_that("conditional chrom", {
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- file.path(tempdir(), "tmp1.gds")
    invisible(file.copy(seqExampleFileName("gds"), gdsfile))
    config <- c(gds_file=file.path(tempdir(), "tmp .gds"),
                phenotype_file=system.file("data", "sample_annotation.RData", package="TopmedPipeline"),
                conditional_variants="1 2",
                conditional_chrom="1",
                n_pcs=0,
                outcome="outcome",
                covars="sex",
                sample_include_file=NA)
    phen <- getPhenotypes(config)
    expect_true(all(c("var_1", "var_2") %in% varLabels(phen$annot)))
    expect_true(all(c("var_1", "var_2") %in% phen$covars))
    unlink(gdsfile)
})
