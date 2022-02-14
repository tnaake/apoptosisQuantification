#' @name getLoadingsOfMarkers
#' 
#' @title Get PCA loadings of marker proteins 
#' 
#' @description
#' The function \code{getLoadingsOfMarkers} retrieves the loadings of the
#' marker proteins (defined by \code{markers}) found in \code{prot}.
#' 
#' \code{getLoadingsOfMarkers} returns a \code{tibble} of the marker 
#' proteins found in \code{prot} containing the loadings for the principal
#' components.
#' 
#' @details 
#' In case there are ambiguous/multiple assignments of proteins (containing 
#' markers) to a loading value (MaxQuant separates these assignments by ";"), 
#' the proteins that are not in \code{markers$protein} are removed.
#' In case there are multiple assignments left of marker proteins, the
#' mean value of the loading values are returned. 
#'
#' @references 
#' Tanzer et al. (2020): Quantitative and Dynamic Catalogs of Proteins 
#' Released during Apoptotic and Necroptotic Cell Death. 
#' \emph{Cell Reports}, 30, 1260-1270.e5. 10.1016/j.celrep.2019.12.079.
#'
#' @param prot matrix, data.frame, or tibble
#' @param markers tibble
#'
#' @return tibble
#'
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column tibble
#'
#' @export
#'
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#' 
#' ## "apoptosis"
#' markers <- readMarkers(type = "apoptosis")
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' se <- readRDS(f)
#' prot <- assay(se) |>
#'     imputeAssay(method = "MinDet")
#' getLoadingsOfMarkers(prot = prot, markers = markers)
#' 
#' ## "nectroptosis"
#' markers <- readMarkers(type = "necroptosis")
#' getLoadingsOfMarkers(prot = prot, markers = markers)
getLoadingsOfMarkers <- function(prot, markers) {
    
    if (!methods::is(markers, "tbl"))
        stop("'markers' is not a tibble")
    
    ## calculate the loadings
    pca <- prot |> 
        as.matrix() |>
        t() |>
        stats::prcomp(center = TRUE, scale = TRUE)
    pca_rotation <- pca$rotation
    
    ## expand markers, some of the features are separated by ";", 
    ## split the strings and replicate the entries
    markers_protein <- splitNames(markers$protein, na.rm = FALSE)#strsplit(markers$protein, split = ";")
    markers_protein_len <- unlist(lapply(markers_protein, length))
    markers <- tibble::tibble(protein = unlist(markers_protein),
        change = rep(markers$change, times = markers_protein_len),
        significant_n = rep(markers$significant_n, times = markers_protein_len))
    
    ## remove duplicated entries
    markers <- markers[!duplicated(markers$protein), ]

    ## look up the names from pca_rotation and calculate the contribution of 
    ## the protein:
    
    ## split the names of the proteins (prot) and remove NA values
    prot_names <- splitNames(protein_names = rownames(pca_rotation), na.rm = TRUE)
    
    ## check if the prot_names are in markers$protein and remove the proteins
    ## that are not in markers
    prot_names <- lapply(prot_names, 
        function(prot_i) prot_i[prot_i %in% markers$protein])
    
    ## get the coordinate in PCn (loadings coordinates)
    marker_vals <- lapply(prot_names, function(prot_names_i) {
        names_i <- prot_names_i[prot_names_i %in% rownames(pca_rotation)]
        if (length(names_i) > 0)
            pca_rotation[names_i, ]
    })
    
    ## take the mean of the product (in case there are several, 
    ## indistinguishable proteins) or return the values as they are
    marker_vals <- lapply(marker_vals, function(vals_i) {
        if (is.matrix(vals_i))
            apply(vals_i, 2, mean)
        else 
            vals_i
    })
    names(marker_vals) <- rownames(pca_rotation)
    
    ## collapse the numeric vectors and create a matrix (the step will 
    ## remove the entries that have a list entry of NULL)
    marker_vals <- do.call("rbind", marker_vals) |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "protein")
    tibble::tibble(marker_vals)
}

#' @name splitNames
#' 
#' @title Split the names when separated by ;
#' 
#' @description
#' The function \code{splitNames} creates a ";"-separated list of protein names. 
#' If \code{na.rm = TRUE}, the \code{NA} values within the protein names 
#' are removed. \code{NA} values might occur if protein names are translated
#' from one id to the other, e.g. from SYMBOL to Ensembl.
#' 
#' @details
#' Helper function for various functions in the \code{apoptosisQuantification}
#' package.
#' 
#' @param protein_names character
#' @param na.rm logical(1)
#' 
#' @examples
#' protein_names <- c("prot_1;prot_2;NA")
#' apoptosisQuantification:::splitNames(protein_names = protein_names,
#'     na.rm = FALSE)
#' apoptosisQuantification:::splitNames(protein_names = protein_names,
#'     na.rm = TRUE)
splitNames <- function(protein_names, na.rm = TRUE) {
    
    ## split the protein names
    protein_names <- strsplit(protein_names, split = ";")
    
    ## remove the NA values
    if (na.rm)
        protein_names <- lapply(protein_names, 
            function(prot_i) prot_i[prot_i != "NA"])
    
    ## return
    protein_names
}


#' @name combineLoadings
#' 
#' @title Combine the loadings/calculate the length of loading vectors
#' 
#' @description
#' The function \code{combineLoadings} calculates the length of loading vectors
#' across several principal components within the Cartesian coordinate system. 
#' 
#' The length of vectors is calculated along the principal coordinates
#' defined by \code{PC}, e.g. if \code{PC = c("PC1", "PC2")} the 
#' length is calculated along the loadings of principal components 1 and 2.  
#' 
#' @details 
#' The length of vectors is calculated using the classical way in the 
#' cartesian coordinate system. 
#' 
#' For 1D (e.g. length of vectors for PC1):
#' \deqn{
#'     \bar{\bar{a}} = \sqrt{a_1^2}
#' }
#' 
#' For 2D (e.g. length of vectors for PC1 and PC2): 
#' \deqn{
#'     \bar{\bar{a}} = \sqrt{a_1^2 + a_2^2}
#' }
#' 
#' For 3D (e.g. length of vectors for PC1, PC2, and PC3):
#' \deqn{
#'     \bar{\bar{a}} = \sqrt{a_1^2 + a_2^2 + a_3^2}
#' }
#' 
#' where $a_1$ is the loading of a data point for PC1, $a_2$ for PC2, $a_3$ for 
#' PC3, etc. 
#' 
#' @param loadings matrix or data.frame
#' @param PC character
#' 
#' @return numeric vector
#' 
#' @importFrom dplyr all_of as_tibble select pull
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble is_tibble
#' 
#' @export
#' 
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#' library(tibble)
#' 
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' se <- readRDS(f)
#' prot <- assay(se)
#' loadings <- prcomp(t(prot))$rotation |> 
#'     as.data.frame() |>
#'     rownames_to_column(var = "protein") |> 
#'     as_tibble()
#' combineLoadings(loadings = loadings, PC = c("PC1", "PC2")) 
combineLoadings <- function(loadings, PC = c("PC1", "PC2")) {
    
    if (!tibble::is_tibble(loadings))
        stop("'loadings' is not a tibble")
    
    ## select the PCs of interest
    names_loadings <- loadings |> 
        dplyr::pull("protein")
    loadings_pcs <- loadings |> 
        tibble::as_tibble() |>
        dplyr::select(dplyr::all_of(PC))
    
    ## square the loadings and calculate the sum
    loadings_pcs <- rowSums(loadings_pcs^2)
    
    ## to obtain the length, take the square root, return the values
    loadings_pcs <- sqrt(loadings_pcs) 
    stats::setNames(loadings_pcs, names_loadings)
}


#' @name createLoadingsTbl
#' 
#' @title Create a tibble containing the loading values for marker and 
#' non-marker proteins
#' 
#' @description
#' The function \code{createLoadingsTbl} creates a tibble containing the loading
#' values for marker and non-marker proteins. 
#' 
#' \code{createLoadingsTbl} returns a tibble with three columns:
#' \itemize{
#'     \item \code{protein}: name of protein (in case of multiple assignemt 
#'         separated by \code{";"}),
#'     \item \code{type}: either \code{"marker"} or \code{non-marker}
#' }
#'
#' @details
#' In case of multiple assignments of proteins, such a multiple protein is
#' treated as a marker protein if there is at least one marker protein
#' (defined in \code{names(vals_markers)}) in the multiple assignment.
#' 
#' @param vals_markers \code{numeric} vector
#' @param vals_all \code{numeric} vector
#' 
#' @return tibble
#' 
#' @export
#' 
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble tibble
#' @importFrom stats setNames
#' 
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#'
#' ## load data set of Tanzer et al. (2020)
#' ## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' 
#' ## get loadings of markers and of the complete data set
#' prot <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis")
#' loadings_markers <- getLoadingsOfMarkers(prot = prot, markers = markers)
#' vals_all <- prcomp(t(prot))$rotation[, 1]
#' vals_markers <- setNames(loadings_markers$PC1, loadings_markers$protein)
#' createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)
createLoadingsTbl <- function(vals_markers = vals_markers, vals_all = vals_all) {
    
    if (!is.vector(vals_markers, mode = "numeric"))
        stop("'vals_markers' is not a numeric vector")
    if (!is.vector(vals_all, mode = "numeric"))
        stop("'vals_all' is not a numeric vector")
    
    ## split the marker proteins
    markers_proteins <- unlist(strsplit(names(vals_markers), split = ";"))

    all_proteins <- splitNames(protein_names = names(vals_all), na.rm = TRUE)
    ind_markers <- unlist(lapply(all_proteins, 
            function (proteins) any(proteins %in% markers_proteins)))

    ind_nonmarkers <- !ind_markers
    vals_nonmarkers <- vals_all[ind_nonmarkers]
    
    ## create data frames
    tbl_markers <- tibble::tibble(
        protein = names(vals_markers), type = "marker", value = vals_markers)
    tbl_nonmarkers <- tibble::tibble(
        protein = names(vals_all)[ind_nonmarkers], type = "non-marker", 
        value = vals_nonmarkers)
    
    ## combine the values and return
    tbl_comb <- rbind(tbl_markers, tbl_nonmarkers)
    tbl_comb
}

#' @name plotECDF
#' 
#' @title Plot empirical cumulative distribution function of loading values
#' 
#' @description 
#' The function \code{plotECDF} creates an empirical cumulative distribution
#' function (ECDF) plot of the loadings distribution of marker and non-markers.
#' 
#' @details 
#' The function \code{plotECDF} takes as input a \code{tibble} (\code{tbl}) 
#' that stores the loading values in the column \code{"value"} and the 
#' information on marker/non-marker in the column \code{"type"}.
#' 
#' @param tbl tibble
#' 
#' @return \code{gg}
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes_string stat_ecdf xlab ylab theme element_blank
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble is_tibble
#' 
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#' 
#' ## load data set of Tanzer et al. (2020)
#' ## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' 
#' ## get loadings of markers and of the complete data set
#' prot <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis")
#' loadings_markers <- getLoadingsOfMarkers(prot = prot, markers = markers)
#' vals_all <- prcomp(t(prot))$rotation[, 1]
#' vals_markers <- setNames(loadings_markers$PC1, loadings_markers$protein)
#' tbl <- createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)
#' 
#' ## plot the distribution as ECDF
#' plotECDF(tbl)
plotECDF <- function(tbl) {
    
    ## check the object
    if (!tibble::is_tibble(tbl)) stop("'tbl' is not a tibble")
    cols <- colnames(tbl)
    if (!"protein" %in% cols) stop("'tbl' does not have column 'protein'")
    if (!"type" %in% cols) stop("'tbl' does not have column 'type'")
    if (!"value" %in% cols) stop("'tbl' does not have column 'value'")
    if (!is.numeric(tbl$value)) stop("'tbl$value' has to be numeric")
    
    ## do the actual plotting
    ggplot2::ggplot(data = tbl, 
        ggplot2::aes_string(x = "value", group = "type", col = "type")) +
        ggplot2::stat_ecdf(size = 1) + 
        ggplot2::xlab("loadings") + ggplot2::ylab("ECDF") + 
        ggplot2::theme(
            legend.title = ggplot2::element_blank(), legend.position = "top")
}

#' @name plotHistogram
#' 
#' @title Plot the histogram (distribution) of loading values 
#' 
#' @description 
#' The function \code{plotHistogram} creates (overlaid) histograms
#' plots of the loadings distribution of marker and non-markers.
#' 
#' @details 
#' The function \code{plotHistogram} takes as input a \code{tibble} (\code{tbl}) 
#' that stores the loading values in the column \code{"value"} and the 
#' information on marker/non-marker in the column \code{"type"}.
#' 
#' @param tbl tibble
#' 
#' @export
#' 
#' @return \code{gg}
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes_string geom_histogram
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay
#' 
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#' 
#' ## load data set of Tanzer et al. (2020)
#' ## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
#' f <- system.file("protein_datasets/tanzer2020.RDS",
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' 
#' ## get loadings of markers and of the complete data set
#' prot <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis", fc = 2, n = 1)
#' loadings_markers <- getLoadingsOfMarkers(prot = prot, markers = markers)
#' vals_all <- prcomp(t(prot))$rotation[, 2]
#' vals_markers <- setNames(loadings_markers$PC2, loadings_markers$protein)
#' tbl <- createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)
#' 
#' plotHistogram(tbl)
plotHistogram <- function(tbl) {
    
    
    ## check the object
    if (!tibble::is_tibble(tbl)) stop("'tbl' is not a tibble")
    cols <- colnames(tbl)
    if (!"protein" %in% cols) stop("'tbl' does not have column 'protein'")
    if (!"type" %in% cols) stop("'tbl' does not have column 'type'")
    if (!"value" %in% cols) stop("'tbl' does not have column 'value'")
    if (!is.numeric(tbl$value)) stop("'tbl$value' has to be numeric")
    
    ## do the actual plotting
    tbl <- tbl |> 
        dplyr::mutate(value = abs(!!as.symbol("value")))
    ggplot2::ggplot(data = tbl, 
        ggplot2::aes_string(x = "value", col = "type", fill = "type")) +
        ggplot2::geom_histogram(alpha = 0.5, position = "identity")
}
    
#' @name performWilcoxonTest
#' 
#' @title Perform the Wilcoxon rank sum and signed rank tests on loading values
#' 
#' @description 
#' The function \code{performWilcoxonTest} performs a Wilcoxon rank sum and 
#' signed rank test on the loading values. The null hypothesis is that 
#' the mean of ranks of the marker loadings is smaller or equal to the 
#' ranks of the non-marker loadings. 
#' The null hypothesis can be rejected with probability \code{p}
#' (alternative hypothesis: loadings of markers are greater than the loadings
#' of non-markers). 
#' 
#' @details 
#' The function \code{performWilcoxonTest} takes as input a \code{tibble} 
#' (\code{tbl}) that stores the loading values in the column \code{"value"} 
#' and the information on marker/non-marker in the column \code{"type"}.
#' 
#' @param tbl tibble
#' 
#' @return htest
#' 
#' @export
#' 
#' @importFrom dplyr filter pull
#' @importFrom MatrixQCvis transformAssay imputeAssay
#' @importFrom stats wilcox.test setNames
#' @importFrom SummarizedExperiment assay
#' 
#' @examples 
#' library(MatrixQCvis)
#' library(SummarizedExperiment)
#' 
#' ## load data set of Tanzer et al. (2020)
#' ## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' 
#' ## get loadings of markers and of the complete data set
#' prot <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis")
#' loadings_markers <- getLoadingsOfMarkers(prot = prot, markers = markers)
#' vals_all <- prcomp(t(prot))$rotation[, 1]
#' vals_markers <- setNames(loadings_markers$PC1, loadings_markers$protein)
#' tbl <- createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)
#' 
#' performWilcoxonTest(tbl)
performWilcoxonTest <- function(tbl) {
    
    ## check the object
    if (!tibble::is_tibble(tbl)) stop("'tbl' is not a tibble")
    cols <- colnames(tbl)
    if (!"protein" %in% cols) stop("'tbl' does not have column 'protein'")
    if (!"type" %in% cols) stop("'tbl' does not have column 'type'")
    if (!"value" %in% cols) stop("'tbl' does not have column 'value'")
    if (!is.numeric(tbl$value)) stop("'tbl$value' has to be numeric")
    
    ## do the actual test
    value_markers <- tbl |>
        dplyr::filter(!!as.symbol("type") == "marker") |>
        dplyr::pull("value")
    value_markers <- abs(value_markers)
    value_nonmarkers <- tbl |>
        dplyr::filter(!!as.symbol("type") == "non-marker") |> 
        dplyr::pull("value")
    value_nonmarkers <- abs(value_nonmarkers)

    ## perform the Wilcoxon test and return the test results
    ## alternative hypothesis is that the value of markers is greater than
    ## the value of non-markers
    wilcox.test(value_markers, value_nonmarkers, alternative = "greater")
}

