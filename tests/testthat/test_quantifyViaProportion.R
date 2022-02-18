set.seed(2022)

## create random protein intensity matrix
prot_vals <- runif(1000, min = 10000, max = 1000000)
prot <- matrix(prot_vals, ncol = 10, nrow = 100)

## fake some protein names and sample names
rownames(prot) <- paste("prot", seq_len(nrow(prot)), sep = "_")
colnames(prot) <- paste("sample", seq_len(ncol(prot)), sep = "_")

inds_mito <- c(28, 37, 87, 69, 29, 12, 27, 13, 54, 52)
prot_mito <- c("MDH2", "OMA1", "NLN", "CS", "OPA1", "CSPG5", "SPATA9", 
     "POLG", "RAB38", "NSD3")
rownames(prot)[inds_mito] <- prot_mito

## function calculateProportionOfMitochondrialProteins
test_that("calculateProportionOfMitochondrialProteins works.", {
    
    prop <- calculateProportionOfMitochondrialProteins(values = prot, 
        id = "Symbol")
    expect_equal(length(prop), 10)
    expect_equal(names(prop), colnames(prot))
    expect_equal(as.vector(prop), 
        c(8.248111, 9.344613, 10.819496, 5.956062, 11.959543, 10.370720,
            7.144944, 10.102232, 11.072084, 9.893684), tolerance = 1e-06)
    prop_foo <- calculateProportionOfMitochondrialProteins(values = prot, 
        id = "Ensembl")
    expect_equal(length(prop), 10)
    expect_equal(names(prop), colnames(prot))
    expect_equal(as.vector(prop_foo), rep(0, 10))
    expect_error(calculateProportionOfMitochondrialProteins(values = 1:10),
        "non-character argument")
    mat_foo <- matrix(1:100, ncol = 10)
    expect_error(calculateProportionOfMitochondrialProteins(values = mat_foo),
        "non-character argument")
    expect_error(calculateProportionOfMitochondrialProteins(values = prot, 
        id = "foo"), "'arg' should be one of ")
})

## function plotProportionOfMitochondrialProteins
test_that("plotProportionOfMitochondrialProteins works.", {
    prop <- calculateProportionOfMitochondrialProteins(values = prot, 
        id = "Symbol")
    gg <- plotProportionOfMitochondrialProteins(proportion = prop)
    expect_is(gg, "gg")
    expect_error(plotProportionOfMitochondrialProteins(proportion = NULL),
        "attempt to set an attribute on NULL")
})

