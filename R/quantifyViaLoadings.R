#' @name getLoadingsOfMarkers
#' 
#' @title Get PCA loadings of marker features 
#' 
#' @description
#' The function \code{getLoadingsOfMarkers} retrieves the loadings of the
#' marker features (defined by \code{markers}) of the \code{values} data set.
#' 
#' \code{getLoadingsOfMarkers} returns a \code{tibble} of the marker 
#' features found in \code{values} containing the loadings for the principal
#' components.
#' 
#' @details 
#' In case there are ambiguous/multiple assignments of features (containing 
#' markers) to a loading value (MaxQuant separates these assignments by ";"), 
#' the features that are not in \code{markers$feature} are removed.
#' In case there are multiple assignments left of marker features, the
#' mean value of the loading values are returned. 
#'
#' @references 
#' Tanzer et al. (2020): Quantitative and Dynamic Catalogs of Proteins 
#' Released during Apoptotic and Necroptotic Cell Death. 
#' \emph{Cell Reports}, 30, 1260-1270.e5. 10.1016/j.celrep.2019.12.079.
#'
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
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
#' apoptosisQuantification:::getLoadingsOfMarkers(values = prot, 
#'     markers = markers)
#'
#' ## "nectroptosis"
#' markers <- readMarkers(type = "necroptosis")
#' apoptosisQuantification:::getLoadingsOfMarkers(values = prot, 
#'     markers = markers)
getLoadingsOfMarkers <- function(values, markers) {
    
    if (!methods::is(markers, "tbl"))
        stop("'markers' is not a tibble")
    
    ## calculate the loadings
    pca <- values |> 
        as.matrix() |>
        t() |>
        stats::prcomp(center = TRUE, scale = TRUE)
    pca_rotation <- pca$rotation
    
    ## expand markers, some of the features are separated by ";", 
    ## split the strings and replicate the entries
    markers_features <- splitNames(markers$feature, na.rm = FALSE)
    markers_features_len <- unlist(lapply(markers_features, length))
    
    markers <- tibble::tibble(feature = unlist(markers_features))
    if ("change" %in% colnames(markers)) {
        markers <- tibble::tibble(markers, 
            change = rep(markers$change, times = markers_features_len))
    }
    if ("significant_n" %in% colnames(markers)) {
        markers <- tibble::tibble(significant_n = rep(markers$significant_n, 
            times = markers_features_len))
    }

    ## remove duplicated entries
    markers <- markers[!duplicated(markers$feature), ]

    ## look up the names from pca_rotation and calculate the contribution of 
    ## the feature:
    
    ## split the names of the features (values) and remove NA values
    features_names <- splitNames(feature_names = rownames(pca_rotation), 
        na.rm = TRUE)
    
    ## check if the feature_names are in markers$feature and remove the features
    ## that are not in markers
    features_names <- lapply(features_names, 
        function(features_i) features_i[features_i %in% markers$feature])
    
    ## get the coordinate in PCn (loadings coordinates)
    marker_vals <- lapply(features_names, function(features_names_i) {
        names_i <- features_names_i[features_names_i %in% rownames(pca_rotation)]
        if (length(names_i) > 0)
            pca_rotation[names_i, ]
    })
    
    ## take the mean of the product (in case there are several, 
    ## indistinguishable features) or return the values as they are
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
        tibble::rownames_to_column(var = "feature")
    tibble::tibble(marker_vals)
}

#' @name splitNames
#' 
#' @title Split the names when separated by ;
#' 
#' @description
#' The function \code{splitNames} creates a ";"-separated list of feature names. 
#' If \code{na.rm = TRUE}, the \code{NA} values within the feature names 
#' are removed. \code{NA} values might occur if feature names are translated
#' from one id to the other, e.g. from SYMBOL to Ensembl.
#' 
#' @details
#' Helper function for various functions in the \code{apoptosisQuantification}
#' package.
#' 
#' @param feature_names character
#' @param na.rm logical(1)
#' 
#' @examples
#' protein_names <- c("prot_1;prot_2;NA")
#' apoptosisQuantification:::splitNames(feature_names = protein_names,
#'     na.rm = FALSE)
#' apoptosisQuantification:::splitNames(feature_names = protein_names,
#'     na.rm = TRUE)
splitNames <- function(feature_names, na.rm = TRUE) {
    
    ## split the feature names
    feature_names <- strsplit(feature_names, split = ";")
    
    ## remove the NA values
    if (na.rm)
        feature_names <- lapply(feature_names, 
            function(names_i) names_i[names_i != "NA"])
    
    ## return
    feature_names
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
#' 
#' loadings <- prcomp(t(prot))$rotation |> 
#'     as.data.frame() |>
#'     rownames_to_column(var = "feature") |> 
#'     as_tibble()
#' 
#' apoptosisQuantification:::combineLoadings(loadings = loadings, 
#'     PC = c("PC1", "PC2")) 
combineLoadings <- function(loadings, PC = c("PC1", "PC2")) {
    
    if (!tibble::is_tibble(loadings))
        stop("'loadings' is not a tibble")
    
    ## select the PCs of interest
    names_loadings <- loadings |> 
        dplyr::pull("feature")
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
#' non-marker features
#' 
#' @description
#' The function \code{createLoadingsTbl} creates a tibble containing the loading
#' values for marker and non-marker features. 
#' 
#' \code{createLoadingsTbl} returns a tibble with three columns:
#' \itemize{
#'     \item \code{feature}: name of feature (e.g. protein, in case of multiple
#'          assignment separated by \code{";"}),
#'     \item \code{type}: either \code{"marker"} or \code{non-marker}
#' }
#'
#' @details
#' In case of multiple assignments of features (e.g. proteins), such a multiple 
#' assignment is treated as a marker feature if there is at least one marker 
#' feature (defined in \code{names(vals_markers)}) in the multiple assignment.
#' 
#' @param vals_markers \code{numeric} vector
#' @param vals_all \code{numeric} vector
#' 
#' @return tibble
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
#' loadings_markers <- apoptosisQuantification:::getLoadingsOfMarkers(
#'     values = prot, markers = markers)
#' vals_all <- prcomp(t(prot))$rotation[, 1]
#' vals_markers <- setNames(loadings_markers$PC1, loadings_markers$feature)
#' apoptosisQuantification:::createLoadingsTbl(vals_markers = vals_markers, 
#'     vals_all = vals_all)
createLoadingsTbl <- function(vals_markers = vals_markers, vals_all = vals_all) {
    
    if (!is.vector(vals_markers, mode = "numeric"))
        stop("'vals_markers' is not a numeric vector")
    if (!is.vector(vals_all, mode = "numeric"))
        stop("'vals_all' is not a numeric vector")
    
    ## split the marker features
    markers_features <- unlist(strsplit(names(vals_markers), split = ";"))

    all_features <- splitNames(feature_names = names(vals_all), na.rm = TRUE)
    ind_markers <- unlist(lapply(all_features, 
            function (features) any(features %in% markers_features)))

    ind_nonmarkers <- !ind_markers
    vals_nonmarkers <- vals_all[ind_nonmarkers]
    
    ## create data frames
    tbl_markers <- tibble::tibble(
        feature = names(vals_markers), type = "marker", value = vals_markers)
    tbl_nonmarkers <- tibble::tibble(
        feature = names(vals_all)[ind_nonmarkers], type = "non-marker", 
        value = vals_nonmarkers)
    
    ## combine the values and return
    tbl_comb <- rbind(tbl_markers, tbl_nonmarkers)
    tbl_comb
}

#' @name prepareTbl
#' 
#' @title Helper function to create the output of createLoadingsTbl
#' 
#' @description
#' The function \code{prepareTbl} prepares the input for 
#' \code{createLoadingsTbl} using the intensity values (\code{values}) and
#' the \code{tibble} with markers (\code{markers}).
#' 
#' @details 
#' The function will call internally the helper function 
#' \code{combineLoadings}. This means the function will return absolute 
#' loadings values even in the case of only one selected principal component.
#' 
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
#' @param PC character, has to be in the format \code{"PC1", "PC2", "PC3", ...},
#' the specification of multiple principal components is possible and the
#' length of loadings vectors will be calculated using these principal 
#' components
#' 
#' @importFrom methods is
#' @importFrom tibble rownames_to_column is_tibble
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
#' 
#' apoptosisQuantification:::prepareTbl(values = prot, markers = markers, 
#'     PC = "PC1")
prepareTbl <- function(values, markers, PC = "PC1") {
    
    if (!(is.matrix(values) | is.data.frame(values)))
        stop("'values' is not a matrix or a data.frame")
    if (!methods::is(markers, "tbl"))
        stop("'markers' is not a tbl")
    
    if (!any(markers$feature %in% rownames(values)))
        stop("There is no 'feature' in 'markers' that was found in 'rownames(values)'")
    
    if (!is.vector(PC, mode = "character"))
        stop("'PC' is not a character vector")
    
    ## create vector with combined loadings of markers for all features
    ## (these will be absolute values)
    loadings <- prcomp(t(values))$rotation |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "feature") |>
        tibble::as_tibble()
    vals_all <- combineLoadings(loadings = loadings, PC = PC) 
    
    ## create vector with combined loadings of markers (these will be absolute 
    ## values)
    loadings_markers <- getLoadingsOfMarkers(values = values, markers = markers)
    vals_markers <- combineLoadings(loadings = loadings_markers, PC = PC)
    
    ## create the tibble by calling the createLoadingsTbl function, this will
    ## create a tibble with the values of markers and non-markers loadings,
    ## the values will be absolute values
    createLoadingsTbl(vals_markers = vals_markers, vals_all = vals_all)
}

#' @name plotPCAandLoadings
#' 
#' @title Plot PCA and loadings of markers and non-markers
#' 
#' @description
#' The function creates a PCA plot and a plot of the loadings. The markers
#' and non-markers are colour-coded in the loadings plot.
#' 
#' @details
#' The function \code{plotPCAandLoadings} takes as input a \code{matrix}, 
#' \code{data.frame}, or \code{tibble} (\code{values}) 
#' that stores the intensities of samples column \code{"value"} and the 
#' information on marker/non-marker in the column \code{"type"}.
#' 
#' @param values \code{matrix} or \code{data.frame}
#' @param markers \code{tibble}
#' @param x_coord \code{character(1)}, specifying the principal component on the
#' x axis
#' @param y_coord \code{character(1)}, specifying the principal component on the
#' y axis
#' 
#' @return \code{plotly}
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab theme_bw
#' @importFrom plotly ggplotly subplot
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
#' ## get the intensities and markers
#' values <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis")
#'
#' ## plot the PCA and loadings plot
#' plotPCAandLoadings(values = values, markers = markers)
plotPCAandLoadings <- function(values, markers, x_coord = "PC1", y_coord = "PC2") {
    
    if (!(is.matrix(values) | is.data.frame(values)))
        stop("'values' is not a matrix or a data.frame")
    if (!tibble::is_tibble(markers))
        stop("'markers' is not a tibble")
    
    if (!any(markers$feature %in% rownames(values)))
        stop("There is no 'feature' in 'markers' that was found in 'rownames(values)'")
    
    if (!is.vector(x_coord, mode = "character") & length(x_coord) != 1)
        stop("'x_coord' is not a character vector of length 1")
    
    if (!is.vector(y_coord, mode = "character") & length(y_coord) != 1)
        stop("'y_coord' is not a character vector of length 1")
    
    loadings_markers <- getLoadingsOfMarkers(values = values, markers = markers)
    
    ## run the PCA
    pca_prcomp <- prcomp(t(values), scale = TRUE, center = TRUE)

    ## get the coordinates
    coords_x <- pca_prcomp$x[, x_coord]
    coords_y <- pca_prcomp$x[, y_coord]

    ## create a data.frame that contains the name of samples (name) and the
    ## x/y coordinates of PCA
    df_coords <- data.frame(name = names(coords_x), values_x = coords_x,
        values_y = coords_y)

    ## loadings markers
    loadings_marker_x <- setNames(loadings_markers[[x_coord]],
        loadings_markers$feature)
    loadings_marker_y <- setNames(loadings_markers[[y_coord]],
        loadings_markers$feature)

    ## get loadings of all features
    loadings_x <- pca_prcomp$rotation[, x_coord]
    loadings_y <- pca_prcomp$rotation[, y_coord]

    ## loadings all
    loadings_x <- createLoadingsTbl(vals_markers = loadings_marker_x,
        vals_all = loadings_x)
    loadings_y <- createLoadingsTbl(vals_markers = loadings_marker_y,
        vals_all = loadings_y)

    ## create data frame that contains the name of features, the type
    ## (marker/non-marker) and the x/y coordinates,
    df_loadings <- data.frame(name = loadings_x[["feature"]],
        type = loadings_x[["type"]], values_x = loadings_x[["value"]],
        values_y = loadings_y[["value"]])

    ## do the plotting of the new coordinates in x_coord and y_coord
    gg_pca <- ggplot2::ggplot(df_coords,
        ggplot2::aes_string(text = deparse(quote(name)))) +
        ggplot2::geom_point(ggplot2::aes_string(x = "values_x", y = "values_y")) +
        ggplot2::xlab(x_coord) + ggplot2::ylab(y_coord) + ggplot2::theme_bw()
    gg_loadings <- ggplot2::ggplot(df_loadings,
        ggplot2::aes_string(text = deparse(quote(name)))) +
        ggplot2::geom_point(ggplot2::aes_string(x = "values_x", y = "values_y"),
            col = "darkgrey", alpha = 0.4,
            data = function(df) subset(df, !!as.symbol("type") == "non-marker")) +
        ggplot2::geom_point(ggplot2::aes_string(x = "values_x", y = "values_y"),
            col = "red", alpha = 0.4,
            data = function(df) subset(df, !!as.symbol("type") == "marker")) +
        ggplot2::xlab(x_coord) + ggplot2::ylab(y_coord) + ggplot2::theme_bw()
    plotly::subplot(
        plotly::ggplotly(gg_pca, tooltip = "text"),
        plotly::ggplotly(gg_loadings, tooltip = "text"),
        titleX = TRUE, titleY = TRUE, nrows = 2, margin = 0.1)
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
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
#' @param PC character, has to be in the format \code{"PC1", "PC2", "PC3", ...},
#' the specification of multiple principal components is possible and the
#' length of loadings vectors will be calculated using these principal 
#' components
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
#' 
#' ## plot the distribution as ECDF
#' plotECDF(values = prot, markers = markers, PC = "PC1")
plotECDF <- function(values, markers, PC = "PC1") {
    
    ## run prepareTbl that will prepare the data and run createLoadingsTbl
    tbl <- prepareTbl(values = values, markers = markers, PC = PC)

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
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
#' @param PC character, has to be in the format \code{"PC1", "PC2", "PC3", ...},
#' the specification of multiple principal components is possible and the
#' length of loadings vectors will be calculated using these principal 
#' components
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
#' ## get intensities and apoptosis markers 
#' values <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis", fc = 2, n = 1)
#' 
#' plotHistogram(values = values, markers = markers, PC = "PC1")
plotHistogram <- function(values, markers, PC = "PC1") {
    
    ## run prepareTbl that will prepare the data and run createLoadingsTbl
    tbl <- prepareTbl(values = values, markers = markers, PC = PC)
    
    ## do the actual plotting
    tbl <- tbl |> 
        dplyr::mutate(value = abs(!!as.symbol("value")))
    ggplot2::ggplot(data = tbl, 
        ggplot2::aes_string(x = "value", col = "type", fill = "type")) +
        ggplot2::geom_histogram(alpha = 0.5, position = "identity")
}

#' @name plotViolin
#' 
#' @title Plot the distribution of loading values via a violin plot
#' 
#' @description 
#' The function \code{plotViolin} creates violin plots of the loadings
#' distribution of marker and non-markers.
#' 
#' @details 
#' The function \code{plotViolin} takes as input a \code{tibble} (\code{tbl}) 
#' that stores the loading values in the column \code{"value"} and the 
#' information on marker/non-marker in the column \code{"type"}.
#' 
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
#' @param PC character, has to be in the format \code{"PC1", "PC2", "PC3", ...},
#' the specification of multiple principal components is possible and the
#' length of loadings vectors will be calculated using these principal 
#' components
#' 
#' @export
#' 
#' @return \code{gg}
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_bw
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
#' values <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis", fc = 2, n = 1)
#' 
#' plotViolin(values = values, markers = markers, PC = "PC1")
plotViolin <- function(values, markers, PC = "PC1") {
    
    ## run prepareTbl that will prepare the data and run createLoadingsTbl
    tbl <- prepareTbl(values = values, markers = markers, PC = PC)

    ## do the actual plotting
    ggplot2::ggplot(data = tbl, 
                    ggplot2::aes_string(y = "value", x = "type")) +
        ggplot2::geom_violin() +
        ggplot2::theme_bw()
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
#' @param values \code{matrix} or \code{data.frame}
#' @param markers tbl
#' @param PC character, has to be in the format \code{"PC1", "PC2", "PC3", ...},
#' the specification of multiple principal components is possible and the
#' length of loadings vectors will be calculated using these principal 
#' components
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
#' values <- assay(tanzer2020)
#' markers <- readMarkers(type = "apoptosis")
#' 
#' performWilcoxonTest(values = values, markers = markers, PC = "PC1")
performWilcoxonTest <- function(values, markers, PC = "PC1") {
    
    ## run prepareTbl that will prepare the data and run createLoadingsTbl
    tbl <- prepareTbl(values = values, markers = markers, PC = PC)
    
    ## prepare the vectors (the first comes from the markers, the second from 
    ## the non-markers)
    value_markers <- tbl |>
        dplyr::filter(!!as.symbol("type") == "marker") |>
        dplyr::pull("value") |>
        abs()
    
    value_nonmarkers <- tbl |>
        dplyr::filter(!!as.symbol("type") == "non-marker") |> 
        dplyr::pull("value") |>
        abs()
    
    ## perform the Wilcoxon test and return the test results
    ## alternative hypothesis is that the value of markers is greater than
    ## the value of non-markers
    wilcox.test(value_markers, value_nonmarkers, alternative = "greater",
        paired = FALSE)
}

