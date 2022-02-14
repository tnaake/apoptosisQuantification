#' @name readMarkers
#' 
#' @title Load a file with known markers linked to apoptosis or necroptosis
#' 
#' @description
#' The function \code{readMarkers} loads a data set with markers linked to
#' apoptosis or necroptosis. 
#' 
#' The returned \code{tibble} contains three columns:
#' \itemize{
#'     \item \code{protein}: name of the marker, 
#'     \item \code{change}: maximum/minimum fold change from the 
#'         time-course experiment (maximum absolute value is taken), 
#'     \item \code{significant_n}: number of significant time-points.
#' }
#' 
#' Classical necrotic markers include lactate deydrogenase (LDH),
#' high-mobility group B1 (HMGB1), myoglobin, enolase, and 14-3-3 proteins.
#' Commonly, a protein release is observed during apoptosis and necroptosis
#' for human myeloid/lymphoma cell lines and human primary macrophage cells as
#' observed time course experiments (cf. Tanzer et al., 2020).  
#' 
#' @details 
#' The function loads the Supplementary Table S1 from Tanzer et al. (2020).
#' 
#' In the data set of Tanzer et al. (2020), the following abbreviations are 
#' used:
#' \itemize{
#'     \item TNF: tumor necrosis factor (TNF-mediated apoptosis),
#'     \item SM: birinapant  (inhibition of cIAPs using small molecules called 
#'         Smac mimetics),
#'     \item IDN: caspase inhibitor IDUN-6556 (IDN-6556, used in induction of 
#'         necroptosis),
#'     \item TNF application: leads to production of cytokines (no induction of 
#'         apoptosis or necroptosis),
#'     \item TNF+SM application: inducton of apptosis,
#'     \item TNF+SM+IDN-6556 application: induction of necroptosis.
#' }
#' 
#' The time points 1 h, 3 h, 5 h, and 7 h are taken into account for defining
#' the markers. For apoptosis, the time points 9 h, 12.5 h and 15 h are 
#' removed prior to filtering the data set. 
#' 
#' @references 
#' Tanzer et al. (2020): Quantitative and Dynamic Catalogs of Proteins 
#' Released during Apoptotic and Necroptotic Cell Death. 
#' \emph{Cell Reports}, 30,  1260-1270.e5. 10.1016/j.celrep.2019.12.079.
#' 
#' @param type character, \code{"apoptosis"} or \code{"necroptosis"}
#' @param fc numeric(1), threshold for fold change (default 2), e.g. if set to
#' 2, the proteins are retained that have absolute fold changes greater or equal
#' to \code{fc}
#' @param n numeric, threshold for number of significant time points
#' (default 1), e.g. if set to 1, the proteins are retained that have at 
#' least 1 significant time point
#' 
#' @return tibble
#' 
#' @importFrom dplyr pull
#' @importFrom tibble tibble
#' 
#' @export
#' 
#' @examples 
#' ## "apoptosis"
#' readMarkers(type = "apoptosis", fc = 2, n = 1)
#' 
#' ## "nectroptosis"
#' readMarkers(type = "necroptosis", fc = 2, n = 2)
readMarkers <- function(type = c("apoptosis", "necroptosis"), fc = 2, n = 1) {
    
    type <- match.arg(type)
    
    ## load the data set that contains the fold changes of the time-courses
    ## sheet "TNF+SM vs co"
    if (type == "apoptosis")
        f <- "protein_markers/Tanzer2020_protein_markers_apoptosis.RDS"
    ## sheet "TNF+SM+IDN-6556 vs co"
    if (type == "necroptosis")
        f <- "protein_markers/Tanzer2020_protein_markers_necroptosis.RDS"
    f <- system.file(f, package = "apoptosisQuantification")
    tbl <- readRDS(f)
    tbl <- tbl[!grepl(x = colnames(tbl), pattern = "9h|12[.]5h|15h")]
    cols <- colnames(tbl)
    
    ## remove 9h, 12.5h, 15h
    cols <- cols[!grepl(x = cols, pattern = "9h|12[.]5h|15h")]

    ## find the columns that contain the information on the significance
    sign_inds <- grep(x = cols, pattern = "Significant.TNF[+]SM")

    ## determine the proteins that are significant in at least n time
    ## course events, depending if type is apoptosis or necroptosis

    sign_logical <- apply(t(tbl[, sign_inds]), 1, 
        function(x) ifelse(is.na(x), FALSE, TRUE))
    sign_n <- apply(sign_logical, 1, sum)

    ## look now into the fold changes and find the absolute highest fold change
    ## for the significant proteins, store these fold changes in a list (fc_max)
    fc_inds <- grep(x = cols, pattern = "Fold.change.[(]Log2[)]")
    fc_tbl <- tbl[, fc_inds]
    fc_max <- lapply(seq_len(nrow(tbl)), function(row_i) {
        fc_tbl_i <- fc_tbl[row_i, ][sign_logical[row_i, ]]
        if (length(fc_tbl_i) > 0) 
            fc_tbl_i[[which.max(abs(fc_tbl_i))]]
        else 
            NA
    })
    fc_max <- unlist(fc_max)

    ## create a logical vector that stores information about which proteins
    ## to take (the ones that have significant times point greater or equal
    ## to n and that have fc_max > 0.5)
    markers_logical <- sign_n >= n & (!is.na(fc_max) & abs(fc_max) >= fc)
    
    ## return the tibble, store inside the protein name, the max (absolute) 
    ## fold change, and the number of time points the fold change was
    ## significant
    tibble::tibble(protein = dplyr::pull(tbl[markers_logical, ], "Gene.names"),
        change = fc_max[markers_logical],
        significant_n = sign_n[markers_logical])
}
