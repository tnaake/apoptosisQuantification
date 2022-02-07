library(MatrixQCvis)
library(SummarizedExperiment)
library(tibble)

markers <- readMarkers(type = "apoptosis", fc = 2, n = 1)
f <- system.file("protein_datasets/tanzer2020.RDS", 
                 package = "apoptosisQuantification")
se <- readRDS(f)
prot <- assay(se)
loadings_markers <- getLoadingsOfMarkers(prot = prot, markers = markers)

## function getLoadingsOfMarkers
test_that("getLoadingsOfMarkers works.", {
    expect_is(loadings_markers, "tbl")
    expect_equal(dim(loadings_markers), c(44, 45))
    expect_equal(colnames(loadings_markers), 
        c("protein", paste("PC", 1:44, sep = "")))
    expect_equal(loadings_markers$protein[1:5], 
        c("SPINT1", "LAMTOR5", "AP1G1", "TIMM8A", "SNRNP200"))
    expect_equal(loadings_markers$PC1[1:5],
        c(0.0060343221, -0.0008232329, -0.0102047172,  0.0008654299, -0.0023743282),
        tolerance = 1e-06)
    expect_equal(loadings_markers$PC2[1:5],
        c(-0.009120662, -0.016252907, -0.035128882, -0.019789178, -0.043396322),
        tolerance = 1e-06)
    
    ## prot
    expect_error(getLoadingsOfMarkers(prot = NULL, markers = markers),
        "'data' must be of a vector type, was 'NULL'")
    expect_error(getLoadingsOfMarkers(prot = c(), markers = markers),
        "'data' must be of a vector type, was 'NULL'")
    mat_foo <- matrix(1:100, ncol = 10)
    expect_error(getLoadingsOfMarkers(prot = mat_foo, 
        markers = markers),
        "non-character argument")
    rownames(mat_foo) <- paste("protein", 1:10, sep = ";")
    expect_equal(dim(getLoadingsOfMarkers(prot = mat_foo, markers = markers)),
        c(0, 1))
    expect_is(getLoadingsOfMarkers(prot = mat_foo, markers = markers), "tbl")
    
    ## markers
    expect_error(getLoadingsOfMarkers(prot = prot, markers = NULL),
        "'markers' is not a tibble")
    expect_error(getLoadingsOfMarkers(prot = prot, markers = c()),
        "'markers' is not a tibble")
    mat_foo <- tibble(protein = c("foo1", "foo2"), fold_change = c(1, 1), 
        significant_n = c(2, 2))
    expect_equal(dim(getLoadingsOfMarkers(prot = prot, markers = mat_foo)),
        c(0, 1))
})

## function splitNames
test_that("splitNames works.", {
    protein_names <- "a;b;c;NA"
    res_1 <- apoptosisQuantification:::splitNames(protein_names = protein_names,
        na.rm = TRUE)
    res_2 <- apoptosisQuantification:::splitNames(protein_names = protein_names,
        na.rm = FALSE)
    expect_is(res_1, "list")
    expect_is(res_2, "list")
    expect_equal(res_1[[1]], c("a", "b", "c"))
    expect_equal(res_2[[1]], c("a", "b", "c", "NA"))
    expect_equal(apoptosisQuantification:::splitNames(protein_names = "a"), 
        list("a"))
    expect_error(apoptosisQuantification:::splitNames(protein_names = 1), 
        "non-character argument")
    expect_error(apoptosisQuantification:::splitNames(protein_names = "a", 
        na.rm = "a"), "argument is not interpretable as logical")
})

## function combineLoadings
loadings <- prcomp(t(prot))$rotation |> 
    as.data.frame() |>
    rownames_to_column(var = "protein") |> 
    as_tibble()
 
test_that("combineLoadings works.", {
    
    loadings_combined <- combineLoadings(loadings = loadings, PC = c("PC1"))
    expect_true(is(loadings_combined, "vector"))
    expect_is(loadings_combined, "numeric")
    expect_equal(length(loadings_combined), 2472)
    expect_equal(length(names(loadings_combined)), 2472)
    
    expect_equal(as.numeric(loadings_combined)[1:5],
        c(0.010278996, 0.013233433, 0.002438168, 0.04417880, 0.013002928),
        tolerance = 1e-06)
    expect_equal(names(loadings_combined)[1:5],
        c("NA;GATD3;NA", "UBA6", "DENND3", "CNOT1", "PGP"))
    expect_equal(sum(loadings_combined), 36.7907, tolerance = 1e-06)
    loadings_combined2 <- combineLoadings(loadings = loadings, 
        PC = c("PC1", "PC2"))
    expect_equal(as.numeric(loadings_combined2)[1:5],
        c(0.024828426, 0.029637857, 0.003028955, 0.052914192, 0.020208540),
        tolerance = 1e-06)
    expect_equal(names(loadings_combined2)[1:5],
        c("NA;GATD3;NA", "UBA6", "DENND3", "CNOT1", "PGP"))
    expect_equal(sum(loadings_combined2), 59.88245, tolerance = 1e-06)
    loadings_combined3 <- combineLoadings(loadings = loadings, 
        PC = c("PC1", "PC2", "PC3"))
    expect_equal(as.numeric(loadings_combined3)[1:5],
        c(0.02561837, 0.03126908, 0.01052189, 0.05559229, 0.02039545),
        tolerance = 1e-06)
    expect_equal(names(loadings_combined3)[1:5],
        c("NA;GATD3;NA", "UBA6", "DENND3", "CNOT1", "PGP"))
    expect_equal(sum(loadings_combined3), 72.80351, tolerance = 1e-06)
    
    ## loadings
    expect_error(combineLoadings(loadings = NULL, PC = c("PC1", "PC2")),
        "'loadings' is not a tibble")
    expect_error(combineLoadings(loadings = tibble(), PC = c("PC1", "PC2")),
        "Column `protein` doesn't exist")
    tbl_foo <- tibble(protein = character(), PC1 = numeric(), PC2 = numeric())
    expect_equal(
        as.numeric(combineLoadings(loadings = tbl_foo, PC = c("PC1", "PC2"))),
        numeric())
    
    ## PC
    expect_error(combineLoadings(loadings = loadings, PC = NULL),
        "must be the same length as the vector")
    expect_error(combineLoadings(loadings = loadings, PC = character()),
        "must be the same length as the vector")
    expect_error(combineLoadings(loadings = loadings, PC = "foo"),
        "Column `foo` doesn't exist")
})

## function createLoadingsTbl
vals_all <- prcomp(t(prot))$rotation[, 1]
vals_markers <- setNames(loadings_markers$PC1, loadings_markers$protein)
tbl <- createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)

test_that("createLoadingsTbl works.", {
    expect_is(tbl, "tbl")
    expect_equal(dim(tbl), c(2472, 3))
    expect_equal(tbl$protein[1:5], 
        c("SPINT1", "LAMTOR5", "AP1G1", "TIMM8A", "SNRNP200"))
    expect_equal(tbl$type[1:5], 
        c("marker", "marker", "marker", "marker", "marker"))
    expect_equal(as.numeric(tbl$value[1:5]), 
        c(0.0060343221, -0.0008232329, -0.0102047172, 0.0008654299, -0.0023743282),
        tolerance = 1e-06)
    expect_equal(as.vector(table(tbl$type)), c(44, 2428))
    expect_equal(sum(tbl$value), 12.92235, tolerance = 1e-06)
    
    ## vals_markers
    expect_error(createLoadingsTbl(vals_markers = NULL, vals_all = vals_all),
        "'vals_markers' is not a numeric vector")
    vals_foo <- setNames(c(0.1, 0.3, 0.1, 0.2, 0.1), paste("p", 1:5, sep = "_"))
    tbl_foo <- createLoadingsTbl(vals_markers = vals_foo, vals_all = vals_all)
    expect_equal(tbl_foo$protein[1:6], 
        c("p_1", "p_2", "p_3", "p_4", "p_5", "NA;GATD3;NA"))
    expect_equal(tbl_foo$type[1:6], 
        c("marker", "marker", "marker", "marker", "marker", "non-marker"))
    expect_equal(as.numeric(tbl_foo$value[1:6]), 
        c(0.1, 0.3, 0.1, 0.2, 0.1, -0.010270), tolerance = 1e-05)
    
    ## vals_all
    expect_error(createLoadingsTbl(vals_markers = vals_markers, 
        vals_all = NULL), "'vals_all' is not a numeric vector")
    tbl_foo <- createLoadingsTbl(vals_markers = vals_markers, 
        vals_all = vals_foo)
    expect_equal(tail(tbl_foo$protein), 
        c("TNPO3", "p_1", "p_2", "p_3", "p_4", "p_5"))
    expect_equal(tail(tbl_foo$type), 
        c("marker", "non-marker", "non-marker", "non-marker", "non-marker", 
            "non-marker"))
    expect_equal(as.numeric(tail(tbl_foo$value)), 
        c(0.02394914, 0.1, 0.3, 0.1, 0.2, 0.1), tolerance = 1e-06)
})

## function plotECDF
test_that("plotECDF works.", {
    expect_is(plotECDF(tbl), "gg")
    expect_error(plotECDF(NULL), "'tbl' is not a tibble")
    expect_error(plotECDF(tibble()), "'tbl' does not have column 'protein'")
    expect_error(plotECDF(tibble(protein = "a")), 
        "'tbl' does not have column 'type'")
    expect_error(plotECDF(tibble(protein = "a", type = "A")), 
        "'tbl' does not have column 'value'")
    expect_error(plotECDF(tibble(protein = "a", type = "A", value = "a")), 
        "'tbl[$]value' has to be numeric")
    expect_is(plotECDF(tibble(protein = "a", type = "A", value = 1)), "gg")
})

## function plotHistogram
test_that("plotHistogram works.", {
    expect_is(plotHistogram(tbl), "gg")
    expect_error(plotHistogram(NULL), "'tbl' is not a tibble")
    expect_error(plotHistogram(tibble()), "'tbl' does not have column 'protein'")
    expect_error(plotHistogram(tibble(protein = "a")), 
                 "'tbl' does not have column 'type'")
    expect_error(plotHistogram(tibble(protein = "a", type = "A")), 
                 "'tbl' does not have column 'value'")
    expect_error(plotHistogram(tibble(protein = "a", type = "A", value = "a")), 
                 "'tbl[$]value' has to be numeric")
    expect_is(plotHistogram(tibble(protein = "a", type = "A", value = 1)), "gg")
})

## function performWilcoxonTest
test_that("performWilcoxonTest works.", {
    wtest <- performWilcoxonTest(tbl)
    expect_is(wtest, "htest")
    expect_equal(wtest$statistic, setNames(45897, "W"))
    expect_equal(wtest$parameter, NULL)
    expect_equal(wtest$p.value, 0.9454851, tolerance = 1e-06)
    expect_equal(wtest$null.value, setNames(0, "location shift"))
    expect_equal(wtest$alternative, "greater")
    expect_equal(wtest$method, 
        "Wilcoxon rank sum test with continuity correction")
    expect_equal(wtest$data.name, "value_markers and value_nonmarkers")
    expect_equal(wtest$p.value, 0.9454851, tolerance = 1e-06) 
    expect_error(performWilcoxonTest(NULL), "'tbl' is not a tibble")
    expect_error(performWilcoxonTest(tibble()), "'tbl' does not have column 'protein'")
    expect_error(performWilcoxonTest(tibble(protein = "a")), 
                 "'tbl' does not have column 'type'")
    expect_error(performWilcoxonTest(tibble(protein = "a", type = "A")), 
                 "'tbl' does not have column 'value'")
    expect_error(performWilcoxonTest(tibble(protein = "a", type = "A", value = "a")), 
                 "'tbl[$]value' has to be numeric")
})
