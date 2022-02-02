#' @title scoreSamples
#'
#' @description 
#' This function performs an enrichment analysis to test whether a given
#' proteomic sample is contaminated with markers from apoptosis or necroptosis.
#' The signatures are taken from Tanzer et al. (2020).
#'
#' @details 
#' The function is insprired by the work of Aurelien Dugourd and his 
#' plasmaContamination package 
#' (https://github.com/saezlab/plasmaContamination/).
#' 
#' @param df data.frame with row names as gene symbols, and then each column
#' should be a sample or the t-values of a differential analysis
#' 
#' @importFrom decoupleR decouple
#' 
#' @export
#' 
#' @return data.frame, containining normalized enrichment (mean) scores for the
#' contamination signature (the higher the score, the higher the potential
#' contamination)
#' 
#' @references
#' Tanzer et al. (2020): Quantitative and Dynamic Catalogs of Proteins 
#' Released during Apoptotic and Necroptotic Cell Death. 
#' \emph{Cell Reports}, 30,  1260-1270.e5. 10.1016/j.celrep.2019.12.079.
#' 
#' @examples
#' library(SummarizedExperiment)
#' 
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' prot <- assay(tanzer2020)
#'
#' ## use the intensities to calculate scores
#' scoreSamples(prot = prot, contamination = "apoptosis", n_perm = 100)
#' 
#' ## use the (modified) loadings vector to calculate scores
#' loadings <- prcomp(t(prot))$rotation |> 
#'     as.data.frame()
#' scoreSamples(prot = loadings, contamination = "apoptosis", n_perm = 100)
#'     
#' 
scoreSamples <- function(prot, contamination = c("apoptosis", "necroptosis"),
    signatures = markers, ...) {

    ## match the contamination argument
    contamination <- match.arg(contamination)

    ## obtain the contamination (apoptoss, necroptosis) markers and prepare the
    ## object for decoupleR
    args_fct <- list(...)
    args_fct[["contamination"]] <- contamination

    ## match the arguments ... with the arguments of readMarkers
    args <- args_fct[names(args_fct) %in% names(formals("readMarkers"))]
    signatures <- do.call("readMarkers", args)
    signatures <- tibble::add_column(signatures, likelihood = 1, 
        set = contamination, mor = signatures$fold_change)

    if (is.vector(prot)) {
        prot <- prot |> 
            as.data.frame()
        colnames(prot) <- "sample"
    }

    ## z-scale
    if ("scale" %in% names(args_fct))
        if (args_fct[["scale"]]) 
            prot <- scale(prot)

    ## fix multiple assignments: if there are multiple assignments write only
    ## one marker protein as prot_names
    prot_names <- splitNames(protein_names = rownames(prot), na.rm = TRUE)
    ind_markers <- lapply(prot_names, 
        function(prot_names_i) prot_names_i %in% signatures$protein)
    prot_names <- lapply(seq_along(prot_names),
        function(i) {
            ind_i <- which(ind_markers[[i]])
            if (length(ind_i) > 0) prot_names[[i]][min(ind_i)]
            else prot_names[[i]][1]
    })
    prot_names <- unlist(prot_names)

    ## convert the matrix to a data.frame
    prot <- prot |> as.data.frame()
    rownames(prot) <- prot_names

    ## create an empty list to store the results
    score_list <- list()

    ## define the arguments for decouple
    if (!("n_perm" %in% names(args_fct))) 
        args_fct[["n_perm"]] <- 1000

    for(i in seq_len(ncol(prot))) {
        sub_prot <- prot |>
            dplyr::select(all_of(i)) |> 
            na.omit()
        sub_prot$decouplerIsGreat <- sub_prot[, 1]

        ## run decoupleR
        means <- decoupleR::run_wmean(mat = as.matrix(sub_prot), 
            network = signatures, .source = "set", .target = "protein", 
            .mor = "mor", .likelihood = "likelihood", 
            times = args_fct[["n_perm"]])

        ## filter the output of decouple, remove some rows and columns and 
        ## store the object in the list
        means <- means |> 
            dplyr::filter(statistic == "norm_wmean") |>
            dplyr::filter(condition != "decouplerIsGreat") |> 
            dplyr::select(-c(statistic, p_value))
        score_list[[i]] <- means
    }

    ## combine the results (norm_wmean scores) and return the object
    score_prot <- Reduce("rbind", score_list)
    colnames(score_prot)[1] <- "set"
    tibble::tibble(score_prot)
}


#' @name plotSampleScores
#' 
#' @title Plot the contamination scores of samples
#' 
#' @description 
#' 
#' @details 
#' 
#' @param scores tibble
#' 
#' @return \code{gg}
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_minimal theme 
#' @importFrom ggplot2 element_text
#' 
#' @examples
#' library(SummarizedExperiment)
#' library(tibble)
#' 
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' prot <- assay(tanzer2020)
#'
#' ## use the intensities to calculate scores
#' scores <- scoreSamples(prot = prot, contamination = "apoptosis", n_perm = 100)
#' scores$condition <- unlist(lapply(
#'     strsplit(x = scores$condition, split = "LFQ[.]intensity[.]"), "[", 2))
#' scores$treatment <- unlist(lapply(
#'     strsplit(scores$condition, split = "_"), "[", 1))
#' plotSampleScores(scores = scores) + 
#'     ggplot2::facet_wrap(~ treatment, scales = "free_x") +
#'     ggplot2::theme(legend.position = "none")
#'     
#' ## use the (modified) loadings vector to calculate scores
#' loadings <- prcomp(t(prot))$rotation |> 
#'     as.data.frame()
#' scores <- scoreSamples(prot = loadings, contamination = "apoptosis", 
#'     n_perm = 100)
#' plotSampleScores(scores = scores) + 
#'     ggplot2::theme(legend.position = "none")
plotSampleScores <- function(scores) {
    
    ## check scores
    if (!"score" %in% colnames(scores)) 
        stop("column 'score' not in 'scores'")
    if (!"condition" %in% colnames(scores)) 
        stop("column 'condition' not in 'scores'")
    if (!"set" %in% colnames(scores)) 
        stop("column 'set' not in 'scores'")
    
    ## factorize the condition column (used for x-axis)
    scores$condition <- factor(scores$condition, levels = scores$condition)
    
    ## do the actual plotting
    ggplot2::ggplot(scores, 
        ggplot2::aes_string(y = "score", x = "condition", 
            group = "set", fill = "set")) +
        ggplot2::geom_bar(position = "dodge", stat = "identity", 
            color="black") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, 
            angle = 90, hjust = 1))
}

