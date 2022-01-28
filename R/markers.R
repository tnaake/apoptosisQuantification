#' @name readMarkers
#' 
#' @description
#' 
#' @details 
#' The function loads the Supplementary Table S1 from Tanzer et al. (2020).
#' 
#' ## classical necrotic markers lactate deydrogenase (LDH), 
#' ## high-mobility group B1 (HMGB1), myoglobin, enolase, and 14-3-3 proteins
#' ## protein releasure during apoptosis and necroptosis (human myeloid/lymphoma cell line and human primary macrphases cells in
#' ## time course experiments)
#' ## TNF: tumor necrosis factor (TNF-mediated apoptosis)
#' ## SM: birinapant  (inhibition of cIAPs using small molecules called 
#' ## Smac mimetics)
#' ## IDN: caspase inhibitor IDUN-6556 (IDN-6556)
#' ## TNF application: lead to production of cytokines
#' ## TNF+SM application: lead to apptosis
#' ## TNF+SM+IDN-6556 application: lead to necroptosis
#' 
#' @references 
#' Tanzer et al. (2020): Quantitative and Dynamic Catalogs of Proteins 
#' Released during Apoptotic and Necroptotic Cell Death. 
#' \textit{Cell Reports}, 30,  1260-1270.e5. 10.1016/j.celrep.2019.12.079.
#' 
#' @param type character, \code{"apoptosis"} or \code{"necroptosis"}
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
#' readMarkers(type = "apoptosis")
#' 
#' ## "nectroptosis"
#' readMarkers(type = "nectroptosis")
readMarkers <- function(type = c("apoptosis", "necroptosis")) {
    
    type <- match.arg(type)
    
    ## load the data set that contains the fold changes of the time-courses
    ## sheet "TNF+SM vs co"
    if (type == "apoptosis")
        f <- "protein_markers/Tanzer2020_protein_markers_apoptosis.RDS"
    ## sheet "TNF+SM+IDN-6556 vs co"
    if (type == "necroptosis")
        f <- "protein_markers/Tanzer2020_protein_markers_necroptosis.RDS")
    f <- system.file(f, package = "apoptosisQuantification")
    tbl <- readRDS(f)
    cols <- colnames(tbl)
    
    ## find the columns that contain the information on the significance
    sign_inds <- grep(x = cols, pattern = "Significant TNF[+]SM")
    
    ## determine the proteins that are significant in at least two/two time
    ## course events, depending if type is apoptosis or necroptosis
    if (type == "apoptosis") n <- 2
    if (type == "necroptosis") n <- 2
    
    sign_logical <- apply(t(tbl[, sign_inds]), 1, 
        function(x) ifelse(is.na(x), FALSE, TRUE))
    sign_n <- apply(sign_logical, 1, sum)
    
    ## look now into the fold changes and find the absolute highest fold change
    ## for the significant proteins, store these fold changes in a list (fc_max)
    fc_inds <- grep(x = cols, pattern = "Fold change [(]Log2[)]")
    fc <- tbl[, fc_inds]
    fc_max <- lapply(seq_len(nrow(tbl)), function(row_i) {
        fc_i <- fc[row_i, ][sign_logical[row_i, ]]
        if (length(fc_i) > 0) 
            fc_i[[which.max(abs(fc_i))]]
        else 
            NA
    })
    fc_max <- unlist(fc_max)

    ## create a logical vector that stores information about which proteins
    ## to take (the ones that have significant times point greater or equal
    ## to n and that have fc_max > 0.5)
    markers_logical <- sign_n >= n & (!is.na(fc_max) & abs(fc_max) > 1)
    
    ## return the tibble, store inside the protein name, the max (absolute) 
    ## fold change, and the number of time points the fold change was
    ## significant
    tibble::tibble(protein = dplyr::pull(tbl[markers_logical, ], "Gene names"),
        fold_change = fc_max[markers_logical],
        significant_n = sign_n[markers_logical])
}
