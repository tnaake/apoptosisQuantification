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
#' The parameter \code{values} has to have (gene) symbols as \code{row.names}.
#' \code{values} can be either a numeric \code{vector}, \code{matrix}, or 
#' \code{data.frame}. It can hold the normalized 
#' 
#' @param values numeric \code{vector}, \code{matrix}, or \code{data.frame}
#' @param contamination \code{character(1)}, either \code{"apoptosis"} or 
#' \code{"necroptosis"}
#' @param ... further arguments passed to \code{readMarkers}
#'  
#' @importFrom stats na.omit
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
#' scoreSamples(values = prot, contamination = "apoptosis", n_perm = 100)
scoreSamples <- function(values, contamination = c("apoptosis", "necroptosis"),
    ...) {

    ## match the contamination argument
    contamination <- match.arg(contamination)

    ## obtain the contamination (apoptoss, necroptosis) markers and prepare the
    ## object for decoupleR
    args_fct <- list(...)
    args_fct[["type"]] <- contamination

    ## match the arguments ... with the arguments of readMarkers
    args <- args_fct[names(args_fct) %in% names(formals("readMarkers"))]
    if (!("signatures" %in% names(args_fct))) {
        args_fct[["signatures"]] <- do.call("readMarkers", args)
    }
    signatures <- tibble::add_column(args_fct[["signatures"]], likelihood = 1, 
        set = contamination, mor = args_fct[["signatures"]]$fold_change)

    if (is.vector(values)) {
        values <- values |> 
            as.data.frame()
        colnames(values) <- "sample"
    }

    ## z-scale
    if ("scale" %in% names(args_fct))
        if (args_fct[["scale"]]) 
            values <- scale(values)

    ## fix multiple assignments: if there are multiple assignments write only
    ## one marker protein as prot_names
    prot_names <- splitNames(protein_names = rownames(values), na.rm = TRUE)
    ind_markers <- lapply(prot_names, 
        function(prot_names_i) prot_names_i %in% signatures$protein)
    prot_names <- lapply(seq_along(prot_names),
        function(i) {
            ind_i <- which(ind_markers[[i]])
            if (length(ind_i) > 0) prot_names[[i]][min(ind_i)]
            else prot_names[[i]][1]
    })
    prot_names <- unlist(prot_names)

    ## convert the matrix to a data.frame, remove duplicated entries
    values <- values |> as.data.frame()
    prot_names_rem <- prot_names[!duplicated(prot_names)]
    values <- values[!duplicated(prot_names), ]
    rownames(values) <- prot_names_rem

    ## create an empty list to store the results
    score_list <- list()

    ## define the arguments for decouple
    if (!("n_perm" %in% names(args_fct))) 
        args_fct[["n_perm"]] <- 1000

    for(i in seq_len(ncol(values))) {
        sub_values <- values |>
            dplyr::select(all_of(i)) |> 
            stats::na.omit()
        sub_values$decouplerIsGreat <- sub_values[, 1]

        ## run decoupleR
        means <- decoupleR::run_wmean(mat = as.matrix(sub_values), 
            network = signatures, .source = "set", .target = "protein", 
            .mor = "mor", .likelihood = "likelihood", 
            times = args_fct[["n_perm"]])

        ## filter the output of decouple, remove some rows and columns and 
        ## store the object in the list
        means <- means |> 
            dplyr::filter(!!as.symbol("statistic") == "norm_wmean") |>
            dplyr::filter(!!as.symbol("condition") != "decouplerIsGreat") |> 
            dplyr::select(-all_of(c("statistic", "p_value")))
        score_list[[i]] <- means
    }

    ## combine the results (norm_wmean scores) and return the object
    score_values <- Reduce("rbind", score_list)
    colnames(score_values)[1] <- "set"
    tibble::tibble(score_values)
}


#' @name plotSampleScores
#' 
#' @title Plot the contamination scores of samples
#' 
#' @description
#' The function \code{plotSampleScores} creates a barplot that visualizes
#' the contamination scores.
#' 
#' @details
#' The function takes the argument \code{scores} that can be created via 
#' \code{scoreSamples}.
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
#' scores <- scoreSamples(values = prot, contamination = "apoptosis",fc = 2, n = 1, n_perm = 100)
#' scores$treatment <- tanzer2020$treatment
#' plotSampleScores(scores = scores) +
#'     ggplot2::facet_wrap(~ treatment, scales = "free_x") +
#'     ggplot2::theme(legend.position = "none")
#' 
#' ## use the intensities to calculate scores
#' scores <- scoreSamples(values = prot, contamination = "necroptosis", n_perm = 100)
#' scores$treatment <- tanzer2020$treatment
#' plotSampleScores(scores = scores) +
#'     ggplot2::facet_wrap(~ treatment, scales = "free_x") +
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

