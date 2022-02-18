library(MatrixQCvis)
library(SummarizedExperiment)
library(tibble)

markers <- readMarkers(type = "apoptosis", fc = 2, n = 1)
f <- system.file("protein_datasets/tanzer2020.RDS", 
                 package = "apoptosisQuantification")
se <- readRDS(f)
prot <- assay(se)
loadings_markers <- getLoadingsOfMarkers(values = prot, markers = markers)

## function getLoadingsOfMarkers
test_that("getLoadingsOfMarkers works.", {
    expect_is(loadings_markers, "tbl")
    expect_equal(dim(loadings_markers), c(44, 45))
    expect_equal(colnames(loadings_markers), 
        c("feature", paste("PC", 1:44, sep = "")))
    expect_equal(loadings_markers$feature[1:5], 
        c("SPINT1", "LAMTOR5", "AP1G1", "TIMM8A", "SNRNP200"))
    expect_equal(loadings_markers$PC1[1:5],
        c(0.0060343221, -0.0008232329, -0.0102047172,  0.0008654299, -0.0023743282),
        tolerance = 1e-06)
    expect_equal(loadings_markers$PC2[1:5],
        c(-0.009120662, -0.016252907, -0.035128882, -0.019789178, -0.043396322),
        tolerance = 1e-06)
    
    ## values
    expect_error(getLoadingsOfMarkers(values = NULL, markers = markers),
        "'values' is not a matrix or a data.frame")
    expect_error(getLoadingsOfMarkers(values = c(), markers = markers),
        "'values' is not a matrix or a data.frame")
    mat_foo <- matrix(1:100, ncol = 10)
    expect_error(getLoadingsOfMarkers(values = mat_foo, 
        markers = markers),
        "non-character argument")
    rownames(mat_foo) <- paste("protein", 1:10, sep = ";")
    expect_equal(dim(getLoadingsOfMarkers(values = mat_foo, markers = markers)),
        c(0, 1))
    expect_is(getLoadingsOfMarkers(values = mat_foo, markers = markers), "tbl")
    
    ## markers
    expect_error(getLoadingsOfMarkers(values = prot, markers = NULL),
        "'markers' is not a tbl")
    expect_error(getLoadingsOfMarkers(values = prot, markers = c()),
        "'markers' is not a tbl")
    mat_foo <- tibble(feature = c("foo1", "foo2"), fold_change = c(1, 1), 
        significant_n = c(2, 2))
    expect_equal(dim(getLoadingsOfMarkers(values = prot, markers = mat_foo)),
        c(0, 1))
})

## function splitNames
test_that("splitNames works.", {
    protein_names <- "a;b;c;NA"
    res_1 <- splitNames(feature_names = protein_names,
        na.rm = TRUE)
    res_2 <- splitNames(feature_names = protein_names,
        na.rm = FALSE)
    expect_is(res_1, "list")
    expect_is(res_2, "list")
    expect_equal(res_1[[1]], c("a", "b", "c"))
    expect_equal(res_2[[1]], c("a", "b", "c", "NA"))
    expect_equal(splitNames(feature_names = "a"), 
        list("a"))
    expect_error(splitNames(feature_names = 1), 
        "non-character argument")
    expect_error(splitNames(feature_names = "a", 
        na.rm = "a"), "argument is not interpretable as logical")
})

## function combineLoadings
loadings <- prcomp(t(prot))$rotation |> 
    as.data.frame() |>
    rownames_to_column(var = "feature") |> 
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
        "'loadings' is not a tbl")
    expect_error(combineLoadings(loadings = tibble(), PC = c("PC1", "PC2")),
        "Column `feature` doesn't exist")
    tbl_foo <- tibble(feature = character(), PC1 = numeric(), PC2 = numeric())
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
vals_markers <- setNames(loadings_markers$PC1, loadings_markers$feature)
tbl <- createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)

test_that("createLoadingsTbl works.", {
    expect_is(tbl, "tbl")
    expect_equal(dim(tbl), c(2472, 3))
    expect_equal(tbl$feature[1:5], 
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
    expect_equal(tbl_foo$feature[1:6], 
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
    expect_equal(tail(tbl_foo$feature), 
        c("TNPO3", "p_1", "p_2", "p_3", "p_4", "p_5"))
    expect_equal(tail(tbl_foo$type), 
        c("marker", "non-marker", "non-marker", "non-marker", "non-marker", 
            "non-marker"))
    expect_equal(as.numeric(tail(tbl_foo$value)), 
        c(0.02394914, 0.1, 0.3, 0.1, 0.2, 0.1), tolerance = 1e-06)
})

test_that("prepareTbl works.", {
    tbl <- prepareTbl(values = prot, markers = markers, PC = "PC1")
    expect_is(tbl, "tbl")
    expect_equal(dim(tbl), c(2472, 3))
    expect_equal(tbl$feature[1:5], 
        c("SPINT1", "LAMTOR5", "AP1G1", "TIMM8A", "SNRNP200"))
    expect_equal(tbl$type[1:5], 
        c("marker", "marker", "marker", "marker", "marker"))
    expect_equal(as.numeric(tbl$value[1:5]), 
        c(0.0060343221, 0.0008232329, 0.0102047172, 0.0008654299, 0.0023743282),
        tolerance = 1e-06)
    expect_equal(as.vector(table(tbl$type)), c(44, 2428))
    expect_equal(sum(tbl$value), 35.87905, tolerance = 1e-06)
    
    ## values
    expect_error(prepareTbl(values = NULL, markers = markers, PC = "PC1"), 
        "'values' is not a matrix or a data.frame")
    mat <- matrix(rnorm(100), nrow = 10)
    expect_error(prepareTbl(values = mat, markers = markers, PC = "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    rownames(mat) <- paste("prot", 1:10)
    expect_error(prepareTbl(values = mat, markers = markers, PC = "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    colnames(mat) <- paste("sample", 1:10)
    expect_error(prepareTbl(values = mat, markers = markers, PC = "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    
    ## markers
    expect_error(prepareTbl(values = prot, markers = NULL, PC = "PC1"), 
        "'markers' is not a tbl")
    
    ## PC
    expect_error(prepareTbl(values = prot, markers = markers, PC = "foo"), 
        "Can't subset columns that don't exist")
    tbl <- prepareTbl(values = prot, markers = markers, PC = c("PC1", "PC2"))
    expect_is(tbl, "tbl")
    expect_equal(dim(tbl), c(2472, 3))
    expect_equal(tbl$feature[1:5], 
        c("SPINT1", "LAMTOR5", "AP1G1", "TIMM8A", "SNRNP200"))
    expect_equal(tbl$type[1:5], 
        c("marker", "marker", "marker", "marker", "marker"))
    expect_equal(as.numeric(tbl$value[1:5]), 
        c(0.01093616, 0.01627374, 0.03658107, 0.01980809, 0.04346123),
        tolerance = 1e-06)
    expect_equal(as.vector(table(tbl$type)), c(44, 2428))
    expect_equal(sum(tbl$value), 59.51404, tolerance = 1e-06)
})

## function plotPCAandLoadings
test_that("plotPCAandLoadings works.", {
    values <- prot
    expect_is(plotPCAandLoadings(values = values, markers = markers, 
        x_coord = "PC1", y_coord = "PC2"), "plotly")
    expect_error(plotPCAandLoadings(values = NULL, markers = markers), 
        "'values' is not a matrix or a data.frame")
    expect_error(plotPCAandLoadings(values = matrix(), markers = markers), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(plotPCAandLoadings(values = matrix(rnorm(100), nrow = 10), 
        markers = markers), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    mat <- matrix(rnorm(100), nrow = 10)
    rownames(mat) <- paste("prot", 1:10)
    colnames(mat) <- paste("sample", 1:10)
    expect_error(plotPCAandLoadings(values = mat, markers = markers), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(plotPCAandLoadings(values = mat, markers = markers, x_coord = "foo"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(plotPCAandLoadings(values = mat, markers = markers, y_coord = "foo"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
})

## function plotECDF
test_that("plotECDF works.", {
    expect_is(plotECDF(values = prot, markers = markers, PC = "PC1"), "gg")
    expect_is(plotECDF(values = prot, markers = markers, PC = c("PC1", "PC2")), "gg")
    expect_error(plotECDF(values = NULL, markers = markers, PC = "PC1"),
        "'values' is not a matrix or a data.frame")
    expect_error(plotECDF(values = matrix(), markers = markers),
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]")
    expect_error(plotECDF(values = prot, markers = NULL),
        "'markers' is not a tbl")
    expect_error(plotECDF(values = prot, markers = markers, PC = "foo"),
        "doesn't exist")
})

## function plotHistogram
test_that("plotHistogram works.", {
    expect_is(plotHistogram(values = prot, markers = markers, PC = "PC1"), "gg")
    expect_error(plotHistogram(values = NULL, markers = markers, PC = "PC1"), 
        "'values' is not a matrix or a data.frame")
    expect_error(plotHistogram(values = matrix(), markers = markers, "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(plotHistogram(values = prot, markers = NULL, PC = "PC1"),
        "'markers' is not a tbl")
    expect_error(plotHistogram(values = prot, markers = markers, PC = 1),
        "'PC' is not a character vector")
    expect_error(plotHistogram(values = prot, markers = markers, PC = "foo"),
        "doesn't exist")
})

## function plotViolin
test_that("plotHistogram works.", {
    expect_is(plotViolin(values = prot, markers = markers, PC = "PC1"), "gg")
    expect_error(plotViolin(values = NULL, markers = markers, PC = "PC1"), 
        "'values' is not a matrix or a data.frame")
    expect_error(plotViolin(values = matrix(), markers = markers, "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(plotViolin(values = prot, markers = NULL, PC = "PC1"),
        "'markers' is not a tbl")
    expect_error(plotViolin(values = prot, markers = markers, PC = 1),
        "'PC' is not a character vector")
    expect_error(plotViolin(values = prot, markers = markers, PC = "foo"),
        "doesn't exist")
})

## function performWilcoxonTest
test_that("performWilcoxonTest works.", {
    wtest <- performWilcoxonTest(values = prot, markers = markers, PC = "PC1")
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
    expect_error(performWilcoxonTest(values = NULL, markers = markers, 
        PC = "PC1"), "'values' is not a matrix or a data.frame")
    expect_error(performWilcoxonTest(values = matrix(), markers = markers, 
        PC = "PC1"), 
        "There is no 'feature' in 'markers' that was found in 'rownames[(]values[)]'")
    expect_error(performWilcoxonTest(values = prot, markers = NULL, PC = "PC1"),
        "'markers' is not a tbl")
    expect_error(performWilcoxonTest(values = prot, markers = markers, PC = 1),
        "'PC' is not a character vector")
    expect_error(performWilcoxonTest(values = prot, markers = markers, 
        PC = "foo"),
        "doesn't exist")
    
})
