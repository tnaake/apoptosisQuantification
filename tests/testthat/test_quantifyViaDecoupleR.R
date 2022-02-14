
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
    expect_equal(scores_apoptosis$value[1:10], 
        c(-4.294760, -5.046478, -5.386511, -4.493642, -4.960600, -4.382777,
            -4.259926, -4.255789, -4.646624, -3.701237), tolerance = 1e-06)
    expect_equal(scores_necroptosis$value[1:10], 
        c(-3.530285, -4.966132, -4.841445, -5.160954, -5.935904, -6.291239,
            -6.135648, -8.223436, -8.004957, -7.784212),
        tolerance = 1e-06)
    expect_equal(sum(scores_apoptosis$value), -117.152, tolerance = 1e-01)
    expect_equal(sum(scores_necroptosis$value), -105.8398, tolerance = 1e-01)
    
    ## fc = 0.1, n = 1, calls decoupleProteins
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            fc = 0.1, n = 1, n_perm = 10)$value), 
        -43.36388, tolerance = 1e-01)
    signatures <- readMarkers("apoptosis", fc = 0.1, n = 1)
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            n_perm = 10, signatures = signatures)$value), 
        -43.36388, tolerance = 1e-01)
    
    ## fc = 3, n = 1, calls decoupleProteins
    signatures <- readMarkers("apoptosis", fc = 3, n = 1)
    expect_equal(sum(scoreSamples(values = prot[, 1:10], 
            contamination = "apoptosis", n_perm = 10, 
            signatures = signatures)$value), 
        -30.16643, tolerance = 1e-01)
    expect_equal(sum(
        scoreSamples(values = prot[, 1:10], contamination = "apoptosis", 
            fc = 3, n = 1, n_perm = 10)$value), 
            -30.16643, tolerance = 1e-01)
    expect_equal(sum(scores_necroptosis$value, fc = 2, n = 2), 
        -111.3441, tolerance = 1e-01)
    expect_error(scoreSamples(values = prot, contamination = "foo"),
        "'arg' should be one of ")
    
    ## fc = 0.1, n = 1, calls permuteProteins
    signatures <- readMarkers("apoptosis", fc = 0.1, n = 1)[, "protein"]
    scores <- scoreSamples(value = prot[, 1:10], contamination = "apoptosis",
        signatures = signatures)
    expect_equal(dim(scores), c(10000, 3))
    expect_equal(scores$value[1:10], 
        c(1.204699, 1.413787, 1.294592, 1.370809, 1.422407, 1.357078, 1.318221,
          1.198024, 1.333419, 1.167783), tolerance = 1e-06)
    expect_equal(sum(scores$value), 10463.81, tolerance = 1e-06)
    expect_equal(as.numeric(table(scores$condition)), rep(1000, 10))
    expect_equal(names(table(scores$condition)), 
        c("DMSO_1h_2", "DMSO_1h_3", "DMSO_1h_4", "DMSO_3h_2", "DMSO_3h_3",
            "DMSO_3h_4", "DMSO_5h_2", "DMSO_5h_3", "DMSO_5h_4", "DMSO_7h_2"))
    expect_equal(unique(scores$set), "apoptosis")
    
    ## fc = 3, n = 1, calls permuteProteins
    signatures <- readMarkers("apoptosis", fc = 3, n = 1)[, "protein"]
    scores <- scoreSamples(values = prot[, 1:10], contamination = "apoptosis",
        signatures = signatures)
    expect_equal(dim(scores), c(10000, 3))
    expect_equal(scores$value[1:10], 
        c(0.7534569, 0.9021550, 0.8396205, 0.7375037, 1.1080970, 1.0890119,
            0.5419135, 0.6787436, 0.9181498, 0.2323270), tolerance = 1e-06)
    expect_equal(sum(scores$value), 14908.81, tolerance = 1e-06)
    expect_equal(as.numeric(table(scores$condition)), rep(1000, 10))
    expect_equal(names(table(scores$condition)), 
        c("DMSO_1h_2", "DMSO_1h_3", "DMSO_1h_4", "DMSO_3h_2", "DMSO_3h_3",
            "DMSO_3h_4", "DMSO_5h_2", "DMSO_5h_3", "DMSO_5h_4", "DMSO_7h_2"))
    expect_equal(unique(scores$set), "apoptosis")
})

## function decoupleProteins
test_that("decoupleProteins works.", {
    args_fct <- list()
    args_fct[["n_perm"]] <- 10
    args_fct[["signatures"]] <- readMarkers("apoptosis", fc = 2, n = 1)
    prot_df <- prot[, 1:10] |>
        as.data.frame()
    scores <- decoupleProteins(values = prot_df, args_fct = args_fct, 
        contamination = "apoptosis")
    expect_equal(dim(scores), c(10, 3))
    expect_equal(unique(scores$set), "apoptosis")
    expect_equal(as.numeric(table(scores$condition)), rep(1, 10))
    expect_equal(names(table(scores$condition)), 
        c("DMSO_1h_2", "DMSO_1h_3", "DMSO_1h_4", "DMSO_3h_2", "DMSO_3h_3",
            "DMSO_3h_4", "DMSO_5h_2", "DMSO_5h_3", "DMSO_5h_4", "DMSO_7h_2"))
    expect_equal(scores$value, 
        c(-3.965659, -4.358913, -4.893849, -4.133633, -4.349551, -4.224200,
            -4.277996, -4.021718, -4.586841, -3.664370), tolerance = 1e-06)
    expect_equal(sum(scores$value), -42.47673, tolerance = 1e-06)
    
    ## entry signatures
    args_foo <- list()
    args_foo[["n_perm"]] <- 10
    expect_error(decoupleProteins(values = prot_df, args_fct = args_foo, 
            contamination = "apoptosis"),
        "'signatures' not in 'names[(]args_fct[)]'")
    
    ## entry n_perm
    args_foo <- list()
    args_foo[["signatures"]] <- readMarkers("apoptosis", fc = 2, n = 1)
    expect_error(decoupleProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' not in 'names[(]args_fct[)]'")
    args_foo[["n_perm"]] <- NULL
    expect_error(decoupleProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' not in 'names[(]args_fct[)]'")
    args_foo[["n_perm"]] <- "abc"
    expect_error(decoupleProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' has to be numeric of length 1")
    args_foo[["n_perm"]] <- c(1, 2)
    expect_error(decoupleProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' has to be numeric of length 1")
})

## function permuteProteins
test_that("permuteProteins works.", {
    args_fct <- list()
    args_fct[["n_perm"]] <- 10
    args_fct[["signatures"]] <- readMarkers("apoptosis", fc = 2, n = 1)
    prot_df <- prot[, 1:10] |>
        as.data.frame()
    scores <- permuteProteins(values = prot_df, args_fct = args_fct, 
        contamination = "apoptosis", seed = 2022)
    expect_equal(dim(scores), c(100, 3))
    expect_equal(unique(scores$set), "apoptosis")
    expect_equal(as.numeric(table(scores$condition)), rep(10, 10))
    expect_equal(names(table(scores$condition)), 
        c("DMSO_1h_2", "DMSO_1h_3", "DMSO_1h_4", "DMSO_3h_2", "DMSO_3h_3",
            "DMSO_3h_4", "DMSO_5h_2", "DMSO_5h_3", "DMSO_5h_4", "DMSO_7h_2"))
    expect_equal(scores$value[1:10], 
        c(1.127443, 1.430533, 1.318550, 1.338809, 1.446879, 1.377358, 1.247451,
            1.229229, 1.466507, 1.054413), tolerance = 1e-06)
    expect_equal(sum(scores$value), 122.0792, tolerance = 1e-06)
    
    ## entry signatures
    args_foo <- list()
    args_foo[["n_perm"]] <- 10
    expect_error(permuteProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"),
        "'signatures' not in 'names[(]args_fct[)]'")
    
    ## entry n_perm
    args_foo <- list()
    args_foo[["signatures"]] <- readMarkers("apoptosis", fc = 2, n = 1)
    expect_error(permuteProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' not in 'names[(]args_fct[)]'")
    args_foo[["n_perm"]] <- NULL
    expect_error(permuteProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' not in 'names[(]args_fct[)]'")
    args_foo[["n_perm"]] <- "abc"
    expect_error(permuteProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' has to be numeric of length 1")
    args_foo[["n_perm"]] <- c(1, 2)
    expect_error(permuteProteins(values = prot_df, args_fct = args_foo, 
        contamination = "apoptosis"), "'n_perm' has to be numeric of length 1")
})

## function plotSampleScores
test_that("plotSampleScores works.", {
    
    ## barplot
    expect_is(gg_apoptosis, "gg")
    expect_is(gg_necroptosis, "gg")
    expect_error(plotSampleScores(scores = 1:10), 
        "column 'value' not in 'scores'")
    df <- data.frame(value = 1:10)
    expect_error(plotSampleScores(scores = df),
        "column 'condition' not in 'scores'")
    df <- data.frame(value = 1:10, condition = paste("a", 1:10))
    expect_error(plotSampleScores(scores = df),
        "column 'set' not in 'scores'")
    df <- data.frame(value = 1:10, condition = paste("a", 1:10), set = "foo")
    expect_is(plotSampleScores(scores = df), "gg")
    expect_equal(gg_apoptosis$labels$x, "samples")
    expect_equal(gg_necroptosis$labels$x, "samples")
    expect_equal(gg_apoptosis$labels$y, "s.d. to empirical null distribution")
    expect_equal(gg_necroptosis$labels$y, "s.d. to empirical null distribution")
    expect_is(gg_apoptosis$layers[[1]]$stat, "StatIdentity")
    expect_is(gg_necroptosis$layers[[1]]$stat, "StatIdentity")
    
    ## boxplot
    signatures <- readMarkers(type = "apoptosis", fc = 2, n = 1)[, "protein"]
    scores <- scoreSamples(values = prot, contamination = "apoptosis", 
        signatures = signatures, n_perm = 10)
    gg <- plotSampleScores(scores = scores)
    expect_is(gg, "gg")
    expect_equal(gg$labels$x, "samples")
    expect_equal(gg$labels$y, "difference to markers")
    expect_is(gg$layers[[1]]$stat, "StatBoxplot")

    
})

