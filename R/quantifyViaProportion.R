#' @name calculateProportionOfMitochondrialProteins
#' 
#' @title Calculate contribution of mitochondrial proteins
#'
#' @description
#' Apoptotic cells show an increase in the proportion of gene expression
#' which originates from the mitochondria, which occurs after the opening of 
#' the MPT pore as apoptosis begins. This gradually results in a decrease in 
#' the total expression of chromosomal genes as the DNA is degraded. 
#' In late apoptosis, there is typically a drastically reduced total count of 
#' genes in a cell, with a large proportion of the gene expression stemming 
#' from mitochondria. 
#'
#' For transcriptomics, typically, at a threshold of where 10-20% of 
#' transcripts are mitochondrial in origin, the total number of transcripts 
#' drops substantially (<20% of what it was prior to apoptosis).
#'
#' The function \code{calculateProportionOfMitochondrialProteins} calculates
#' proportions of the intensities of mitochondrial proteins (i.e. proteins that 
#' are encoded in the mitochondrial DNA) to all intensities. The function 
#' returns the proportions in percent. 
#'
#' Mitochondrial proteins are defined by the Human Protein Atlas 
#' (1139 proteins found by searching the Human Protein Atlas database for 
#' mitochondrial sublocalization). These proteins are matched against the
#' proteins found in \code{prot}. Currently, two types of ids are 
#' encoded: either \code{rownames(prot)} have to be "Symbol" or "Ensembl" ids.
#'
#' @details 
#' Protein intensities have to in untransformed format such that 
#' proportions can be calculated as are.
#'
#' Protein ids are stored in \code{rownames(prot)}. In case of multiple
#' assignments the ids are separeted by \code{";"}. In that case, a 
#' feature will be treated as of mitochondrial origin when there is at least
#' one mitochondrial protein in the multiple assignment. 
#'
#' @references 
#' The Human Protein Atlas 
#' (search: subcell_location:Mitochondria, 
#' https://www.proteinatlas.org/humanproteome/subcellular/mitochondria, 
#' accessed January 31, 2022).
#' 
#' Uhlen et al. (2015): Tissue-based map of the human proteome.
#' \emph{Science}, 347, 1260419. 10.1126/science.1260419
#' 
#' Thul et al. (2017): A subcellular map of the human proteome.
#' \emph{Science}, 356, eaal3321. 10.1126/science.aal3321 
#'
#' @param prot matrix, data.frame, or tibble
#' @param id character(1), either \code{"Symbol"} or \code{"Ensembl"}
#'
#' @return numeric vector
#'
#' @export
#' 
#' @importFrom utils read.table
#'
#' @examples 
#' set.seed(2022) 
#' ## fake measured protein intensities and create matrix that stores 
#' ## intensities
#' prot_vals <- runif(1000, min = 10000, max = 1000000)
#' prot <- matrix(prot_vals, ncol = 10, nrow = 100)
#' 
#' ## fake some protein names and sample names
#' rownames(prot) <- paste("prot", seq_len(nrow(prot)), sep = "_")
#' colnames(prot) <- paste("sample", seq_len(ncol(prot)), sep = "_")
#' 
#' ## randomly assign 10 mitochondrial proteins to the matrix
#' inds_mito <- sample(seq_len(nrow(prot)), 10)
#' prot_mito <- c("MDH2", "OMA1", "NLN", "CS", "OPA1", "CSPG5", "SPATA9", 
#'     "POLG", "RAB38", "NSD3")
#' rownames(prot)[inds_mito] <- prot_mito
#' 
#' calculateProportionOfMitochondrialProteins(prot = prot, id = "Symbol")
calculateProportionOfMitochondrialProteins <- function(prot, 
    id = c("Symbol", "Ensembl")) {
    
    id <- match.arg(id)
    
    f <- "subcellular_location/subcell_location_Mitochondria.tsv"
    f <- system.file(f, package = "apoptosisQuantification")
    mito_proteins <- utils::read.table(f, sep = "\t", header = TRUE)
    
    if (id == "Symbol") col_id <- "Gene"
    if (id == "Ensembl") col_id <- "Ensembl"
    
    ## get the ids of the mitochondrial proteins of the database
    protein_mito_id <- mito_proteins[, col_id]
    
    ## split the names of prot (separated by ";"), remove NAs
    prot_names <- strsplit(rownames(prot), split = ";")
    prot_names <- lapply(prot_names, 
        function(prot_names_i) prot_names_i[prot_names_i != "NA"])

    ## intersect with the ids from prot
    ## (in case of multiple assignments treat the feature as of mitochondrial
    ## origin if there is at least one mitochondrial protein)
    protein_mito <- lapply(prot_names, 
        function(prot_names_i) any(prot_names_i %in% protein_mito_id))
    protein_mito <- unlist(protein_mito)
    
    ## calculate the sum of intensities for mitochondrial and all proteins
    sum_mito <- prot[protein_mito, , drop = FALSE]
    sum_mito <- apply(sum_mito, 2, sum, na.rm = FALSE)
    sum_all <- apply(prot, 2, sum, na.rm = FALSE)
    
    ## calculate the proportion and return
    setNames(sum_mito / sum_all * 100, colnames(prot))
}

#' @name plotProportionOfMitochondrialProteins
#' 
#' @title Plot proportion of mitochondrial proteins
#' 
#' @description
#' The function \code{plotProportionOfMitochondrialProteins} creates a 
#' plot of the proportion of mitochondrial proteins 
#' (obtained via \code{calculateProportionOfMitochondrialProtein}).
#' 
#' @details 
#' \code{names(proportion)} are taken as x-axis labels. 
#' 
#' If there are no names for \code{proportion} 
#' (\code{names(proportion)} is NULL), names from \code{1:length(proportion)}
#' will be assigned to \code{proportion}.
#'
#' @param proportion numeric vector
#'
#' @return \code{gg}
#'
#' @export
#' 
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw ylab theme 
#' @importFrom ggplot2 element_text
#' 
#' @examples 
#' set.seed(1)
#' ## fake measured protein intensities and create matrix that stores 
#' ## intensities
#' prot_vals <- runif(1000, min = 10000, max = 1000000)
#' prot <- matrix(prot_vals, ncol = 10, nrow = 100)
#' 
#' ## fake some protein names and sample names
#' rownames(prot) <- paste("prot", seq_len(nrow(prot)), sep = "_")
#' colnames(prot) <- paste("sample", seq_len(ncol(prot)), sep = "_") 
#' 
#' ## randomly assign 10 mitochondrial proteins to the matrix
#' inds_mito <- sample(seq_len(nrow(prot)), 10)
#' prot_mito <- c("MDH2", "OMA1", "NLN", "CS", "OPA1", "CSPG5", "SPATA9", "POLG", "RAB38", "NSD3")
#' rownames(prot)[inds_mito] <- prot_mito
#' 
#' prop_mito <- calculateProportionOfMitochondrialProteins(prot = prot, id = "Symbol")
#' plotProportionOfMitochondrialProteins(proportion = prop_mito)
plotProportionOfMitochondrialProteins <- function(proportion) {
    
    ## add 1:length(proportion) as names for proportion if there are no
    ## pre-defined names
    if (is.null(names(proportion)))
        names(proportion) <- seq_len(length(proportion))
    
    proportion <- tibble::tibble(value = proportion, 
        sample = names(proportion))
    proportion$sample <- factor(proportion$sample, levels = proportion$sample)
    
    ## plot the proportions
    ggplot2::ggplot(proportion,
        ggplot2::aes_string(x = "sample", y = "value")) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::ylab("proportion mitochondrial proteins (%)") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}
