
library(SummarizedExperiment)
library(tibble)
 
f <- system.file("protein_datasets/tanzer2020.RDS", 
    package = "apoptosisQuantification")
tanzer2020 <- readRDS(f)
prot <- assay(tanzer2020)

## use the intensities to calculate scores
scores_apoptosis <- scoreSamples(values = prot, contamination = "apoptosis",
    fc = 2, n = 1, n_perm = 10)
scores_apoptosis$treatment <- tanzer2020$treatment
gg_apoptosis <- plotSampleScores(scores = scores_apoptosis) +
    ggplot2::facet_wrap(~ treatment, scales = "free_x") +
    ggplot2::theme(legend.position = "none")

## use the intensities to calculate scores
scores_necroptosis <- scoreSamples(values = prot, contamination = "necroptosis", 
    fc = 2, n = 2, n_perm = 10)
scores_necroptosis$treatment <- tanzer2020$treatment
gg_necroptosis <- plotSampleScores(scores = scores_necroptosis) +
    ggplot2::facet_wrap(~ treatment, scales = "free_x") +
    ggplot2::theme(legend.position = "none")

## function scoreSamples
test_that("scoreSamples works.", {
    expect_is(scores_apoptosis, "tbl")
    expect_equal(dim(scores_apoptosis), c(44, 4))
    expect_equal(as.numeric(table(scores_apoptosis$set)), 44)
    expect_equal(names(table(scores_apoptosis$set)), "apoptosis")
    expect_equal(dim(scores_necroptosis), c(44, 4))
    expect_equal(as.numeric(table(scores_necroptosis$set)), 44)
    expect_equal(names(table(scores_necroptosis$set)), "necroptosis")
    expect_true(is.character(scores_apoptosis$condition))
    expect_true(is.character(scores_necroptosis$condition))
    expect_equal(sum(scores_apoptosis$score), -117.152, tolerance = 1e-01)
    expect_equal(sum(scores_necroptosis$score), -105.8398, tolerance = 1e-01)
    
    ## fc = 0.1, n = 1
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            fc = 0.1, n = 1, n_perm = 10)$score), 
        -43.36388, tolerance = 1e-01)
    signatures <- readMarkers("apoptosis", fc = 0.1, n = 1)
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            n_perm = 10, signatures = signatures)$score), 
        -43.36388, tolerance = 1e-01)
    
    ## fc = 3, n = 1
    signatures <- readMarkers("apoptosis", fc = 3, n = 1)
    expect_equal(sum(scoreSamples(values = prot[, 1:10], 
            contamination = "apoptosis", n_perm = 10, 
            signatures = signatures)$score), 
        -30.16643, tolerance = 1e-01)
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            fc = 3, n = 1, n_perm = 10)$score), 
            -30.16643, tolerance = 1e-01)
    expect_equal(sum(scores_necroptosis$score, fc = 2, n = 2), 
        -111.3441, tolerance = 1e-01)
    expect_error(scoreSamples(values = prot, contamination = "foo"),
        "'arg' should be one of ")
})

## function plotSampleScores
test_that("plotSampleScores works.", {
    expect_is(gg_apoptosis, "gg")
    expect_is(gg_necroptosis, "gg")
    expect_error(plotSampleScores(scores = 1:10), 
        "column 'score' not in 'scores'")
    df <- data.frame(score = 1:10)
    expect_error(plotSampleScores(scores = df),
        "column 'condition' not in 'scores'")
    df <- data.frame(score = 1:10, condition = paste("a", 1:10))
    expect_error(plotSampleScores(scores = df),
        "column 'set' not in 'scores'")
    df <- data.frame(score = 1:10, condition = paste("a", 1:10), set = "foo")
    expect_is(plotSampleScores(scores = df), "gg")
})

