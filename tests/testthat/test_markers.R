## function readMarkers
test_that("readMarkers works.", {
    apoptosis <- readMarkers(type = "apoptosis", fc = 1, n = 1)
    necroptosis <- readMarkers(type = "necroptosis", fc = 1, n = 2)
    expect_is(apoptosis, "tbl")
    expect_equal(colnames(apoptosis), 
        c("protein", "change", "significant_n"))
    expect_equal(dim(apoptosis), c(89, 3))
    expect_equal(apoptosis$protein[1:5], 
        c("ACSL4", "ADNP", "AP1G1", "ASPSCR1", "ASRGL1"))
    expect_equal(apoptosis$change[1:5], 
        c(2.477711, 2.057616, 2.040779, 1.479240, 2.152555),
        tolerance = 1e-06)
    expect_equal(apoptosis$significant_n[1:5], 
        c(1, 1, 1, 1, 1))
    expect_is(necroptosis, "tbl")
    expect_equal(colnames(necroptosis), 
        c("protein", "change", "significant_n"))
    expect_equal(dim(necroptosis), c(253, 3))
    expect_equal(necroptosis$protein[1:5], 
        c("ACIN1", "ACOT13", "ACSL4", "ACY1;ABHD14A-ACY1", "ADAR"))
    expect_equal(necroptosis$change[1:5], 
        c(2.117764, -2.305560, 1.440479, -2.364775, 3.851311),
        tolerance = 1e-06)
    expect_equal(necroptosis$significant_n[1:5], 
        c(3, 2, 2, 2, 3))
    
    ## type
    expect_error(readMarkers(type = "foo", fc = 1, n = 2),
        "'arg' should be one of ")
    
    ## fc
    expect_error(readMarkers(type = "apoptosis", fc = NULL, n = 2),
        "Must subset rows with a valid subscript vector")
    expect_equal(dim(readMarkers(type = "apoptosis", fc = "foo", n = 2)),
        c(0, 3))
    expect_equal(dim(readMarkers(type = "apoptosis", fc = 3, n = 1)), c(20, 3))
    
    ## n
    expect_error(readMarkers(type = "apoptosis", fc = 1, n = NULL),
        "Must subset rows with a valid subscript vector")
    expect_equal(dim(readMarkers(type = "apoptosis", fc = 1, n = "foo")),
        c(0, 3))
    expect_equal(dim(readMarkers(type = "apoptosis", fc = 1, n = 2)), c(5, 3))
})

