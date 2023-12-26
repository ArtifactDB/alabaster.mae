# This tests the stageMultiAssayExperiment function.
# library(testthat); library(alabaster.mae); source("test-MultiAssayExperiment.R")

set.seed(100)
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

library(SummarizedExperiment)
se <- SummarizedExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

library(MultiAssayExperiment)
mae <- MultiAssayExperiment(list(rnaseq=se))

test_that("saveMultiAssayExperiment works as expected", {
    temp <- tempfile()
    saveObject(mae, temp)

    out <- readMultiAssayExperiment(temp)
    expect_identical(sampleMap(mae), sampleMap(out))
    expect_identical(colData(mae), colData(out))
    expect_identical(names(experiments(out)), names(experiments(out)))
})

test_that("stageMultiAssayExperiment works as expected with non-empty colData", {
    colData(mae)$foo <- runif(10)
    temp <- tempfile()
    saveObject(mae, temp)
    out <- readMultiAssayExperiment(temp)
    expect_equal(colData(mae), colData(out))
})

test_that("stageMultiAssayExperiment works as expected with non-trivial sample mappings", {
    mae <- MultiAssayExperiment(list(rnaseq=se[,1:5], chipseq=se[1:20,3:10]))
    temp <- tempfile()
    saveObject(mae, temp)

    out <- readMultiAssayExperiment(temp)
    expect_identical(names(experiments(out)), c("rnaseq", "chipseq"))

    s1 <- sampleMap(mae)
    s1$assay <- as.character(s1$assay)
    s2 <- sampleMap(out)
    s2$assay <- as.character(s2$assay)
    expect_identical(s1, s2)

    # Even more complicated, this time with experiment-specific names.
    exps <- list(rnaseq=se[,1:5], chipseq=se[1:20,3:10], facs=se[5:10, c(1, 2, 9)])
    sample_map <- list()
    for (n in names(exps)) {
        sample.id <- colnames(exps[[n]])
        colnames(exps[[n]]) <- paste0("rna:", sample.id)
        sample_map[[n]] <- DataFrame(primary=sample.id, assay=n, colname=colnames(exps[[n]]))
    }
    sample_map <- do.call(rbind, sample_map)
    sample_map$assay <- factor(sample_map$assay)

    mae <- MultiAssayExperiment(
        exps, 
        colData=DataFrame(row.names=colnames(se), FOO=runif(ncol(se))),
        sampleMap=sample_map
    )
    temp <- tempfile()
    saveObject(mae, temp)

    out <- readMultiAssayExperiment(temp)
    expect_identical(names(experiments(out)), c("rnaseq", "chipseq", "facs"))

    s1 <- sampleMap(mae)
    s1$assay <- as.character(s1$assay)
    s2 <- sampleMap(out)
    s2$assay <- as.character(s2$assay)
    expect_identical(s1, s2)
})

test_that("stageMultiAssayExperiment respects entries in the metadata", {
    metadata(mae)$WHEE <- "foo"
    temp <- tempfile()
    saveObject(mae, temp)
    mae2 <- readMultiAssayExperiment(temp)
    expect_identical(metadata(mae2)$WHEE, "foo")
})
