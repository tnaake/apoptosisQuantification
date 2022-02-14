#' @title scoreSamples
#'
#' @description 
#' This function performs an enrichment analysis to test whether a given
#' proteomic sample is contaminated with markers from apoptosis or necroptosis.
#'
#' @details 
#' The parameter \code{values} has to have (gene) symbols as \code{row.names}.
#' \code{values} can be either a numeric \code{vector}, \code{matrix}, or 
#' \code{data.frame}. It holds the normalized intensity values of proteins.
#' 
#' The signatures are taken per default from Tanzer et al. (2020), but also 
#' a \code{tibble} containing the contamination markers can be given to the 
#' function as argument \code{signatures}. 
#' In that case the \code{tibble} has to contain at least the 
#' column \code{"protein"} that contains the marker proteins.  
#' 
#' The function will either run \code{decoupleProteins} or 
#' \code{permuteProteins} to calculate scores. The behaviour will be triggered
#' by the presence of the column \code{value} in \code{signatures} (in that 
#' case, \code{decoupleProteins} is run).
#' 
#' @param values numeric \code{vector}, \code{matrix}, or \code{data.frame}
#' @param contamination \code{character(1)}, either \code{"apoptosis"} or 
#' \code{"necroptosis"}
#' @param ... further arguments passed to \code{readMarkers}
#'  
#' @importFrom methods formalArgs
#' @importFrom stats na.omit
#' @importFrom decoupleR decouple
#' 
#' @export
#' 
#' @return data.frame, containining normalized enrichment (mean) scores for the
#' contamination signature (the higher the score, the higher the potential
#' contamination). If there is no quantitative information on changes 
#' (e.g. fold changes or t-values), the function will return an empirical
#' distribution of difference of means of the markers to all proteins in the 
#' data set
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
#' ## use the intensities to calculate scores, the function will load the 
#' ## data set of Tanzer et al. (2020)
#' scoreSamples(values = prot, contamination = "apoptosis", n_perm = 100)
#' 
#' ## use the contamination scores of Tanzer et al. (2020), but with specified
#' ## fc and n_significan
scoreSamples <- function(values, contamination = c("apoptosis", "necroptosis"),
    ...) {

    ## match the contamination argument
    contamination <- match.arg(contamination)

    ## obtain the contamination (apoptoss, necroptosis) markers and prepare the
    ## object for decoupleR
    args_fct <- list(...)
    args_fct[["type"]] <- contamination

    ## match the arguments ... with the arguments of readMarkers
    args <- args_fct[names(args_fct) %in% methods::formalArgs("readMarkers")]
    if (!("signatures" %in% names(args_fct))) {
        args_fct[["signatures"]] <- do.call("readMarkers", args)
    }

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
        function(prot_names_i) prot_names_i %in% args_fct[["signatures"]]$protein)
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

    ## define the arguments for decouple
    if (!("n_perm" %in% names(args_fct))) 
        args_fct[["n_perm"]] <- 1000

    ## depening if the column value is present in signatures, either run the
    ## function decoupleProteins or permuteProteins (otherwise)
    if ("change" %in% colnames(args_fct[["signatures"]])) {
        score_values <- decoupleProteins(values = values, args_fct = args_fct,
            contamination = contamination)
    } else {
        score_values <- permuteProteins(values = values, args_fct = args_fct,
            contamination = contamination)
    }

    score_values
}

#' @name decoupleProteins
#' 
#' @title Calculate score values using decoupleR/run_wmean
#' 
#' @description 
#' The function \code{decoupleProteins} calculates contamination scores
#' using the \code{run_wmean} function from the \code{decoupleR}. The 
#' function will return score values (per sample) that can be interpreted as the 
#' standard deviations away from an empirical null distribution.
#' 
#' @details 
#' The function is insprired by the work of Aurelien Dugourd and his 
#' plasmaContamination package 
#' (https://github.com/saezlab/plasmaContamination/).
#'
#' The function will be called in \code{scoreSamples}.
#' 
#' @param values matrix, the rownames contain the protein names (SYMBOL ids) and
#' colnames the samples
#' @param args_fct list, containing the entries \code{"n_perm"} 
#' (\code{numeric(1)} and \code{"signatures"} (\code{tibble})
#' @param contamination \code{character(1)}, either \code{"apoptosis"} or
#' \code{"necroptosis"}
#' 
#' @return tibble
#' 
#' @export
#' 
#' @examples 
#' library(SummarizedExperiment)
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' prot <- assay(tanzer2020) |>
#'     as.data.frame()
#' 
#' ## define a list with paramters
#' args_fct <- list()
#' args_fct[["signatures"]] <- readMarkers(type = "apoptosis", fc = 2, n = 1)
#' args_fct[["n_perm"]] <- 100
#' 
#' ## run the function
#' decoupleProteins(values = prot, args_fct = args_fct, 
#'     contamination = "apoptosis")
decoupleProteins <- function(values, args_fct, 
    contamination = c("apoptosis", "necroptosis")) {
    
    contamination <- match.arg(contamination)
    
    if (!("n_perm" %in% names(args_fct)))
        stop("'n_perm' not in 'names(args_fct)'")
    
    n_perm <- args_fct[["n_perm"]]
    if (!(is.numeric(n_perm) & length(n_perm) == 1))
        stop("'n_perm' has to be numeric of length 1")
    
    if (!("signatures" %in% names(args_fct)))
        stop("'signatures' not in 'names(args_fct)'")
    signatures <- args_fct[["signatures"]]
    
    signatures <- tibble::add_column(signatures, 
        likelihood = 1, set = contamination, mor = signatures$change)
    
    ## create an empty list to store the results
    score_list <- list()
    
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
    score_tbl <- Reduce("rbind", score_list)
    colnames(score_tbl)[colnames(score_tbl) == "source"] <- "set"
    colnames(score_tbl)[colnames(score_tbl) == "score"] <- "value"
    tibble::tibble(score_tbl)
}

#' @name permuteProteins
#' 
#' @title Calculate score values by calculating difference between
#' means of permuted proteins and markers
#' 
#' @description 
#' The function \code{permuteProteins} calculates the difference between the
#' mean of marker proteins and random samples of the protein universe 
#' (proteins in \code{values}). The scores will be a relative quantification
#' of the marker proteins to the universe protein intensities. Samples that 
#' undergo apoptosis (or any other contamination) will have different 
#' differences than uncontaminated samples.
#' 
#' The output will be a \code{tibble} containing the differences in means
#' for the samples and repetitions. 
#' 
#' @details 
#' The number of repetitions is defined by \code{n_perm}.
#' 
#' The function will be called in \code{scoreSamples}.
#' 
#' @param values \code{data.frame}
#' @param args_fct list of arguments, has to contain the entries
#' \code{n_perm} (\code{numeric(1)}) and \code{signatures} (\code{tibble})
#' @param contamination \code{character(1)}, either \code{"apoptosis"} or 
#' \code{"necroptosis"}
#' @param seed \code{numeric(1)}
#'
#' @importFrom dplyr everything
#' @importFrom tidyr pivot_longer
#' 
#' @export
#' 
#' @examples
#' library(SummarizedExperiment)
#' f <- system.file("protein_datasets/tanzer2020.RDS", 
#'     package = "apoptosisQuantification")
#' tanzer2020 <- readRDS(f)
#' prot <- assay(tanzer2020) |>
#'     as.data.frame()
#' 
#' ## define a list with paramters
#' args_fct <- list()
#' args_fct[["signatures"]] <- readMarkers(type = "apoptosis", fc = 2, n = 1)
#' args_fct[["n_perm"]] <- 100
#' 
#' ## run the function
#' permuteProteins(values = prot, args_fct = args_fct, 
#'     contamination = "apoptosis", seed = 2022)
permuteProteins <- function(values, args_fct, 
    contamination = c("apoptosis", "necroptosis"), seed = 2022) {

    contamination <- match.arg(contamination)

    if (!("n_perm" %in% names(args_fct)))
        stop("'n_perm' not in 'names(args_fct)'")

    n_perm <- args_fct[["n_perm"]]
    if (!(is.numeric(n_perm) & length(n_perm) == 1))
        stop("'n_perm' has to be numeric of length 1")

    if (!("signatures" %in% names(args_fct)))
        stop("'signatures' not in 'names(args_fct)'")
    signatures <- args_fct[["signatures"]]

    ## calculate the means of the markers
    mean_marker <- apply(values[signatures$protein, ], 2, mean, na.rm = TRUE)

    ## create random numbers and calculate the means for these random draws
    set.seed(seed)
    inds <- lapply(seq_len(n_perm), function(i) {
        sample(x = seq_len(nrow(values)), size = nrow(signatures), 
               replace = FALSE)
    })
    means_inds <- lapply(inds, function(inds_i) {
          apply(values[inds_i, ], 2, mean, na.rm = TRUE)
    })

    ## calculate the difference between the markers and the randomly sampled
    ## proteins
    diff_means_inds <- lapply(means_inds, function(means_inds_i) {
        means_inds_i - mean_marker
    })

    score_values <- do.call("rbind", diff_means_inds) |>
        as.data.frame() |>
        tidyr::pivot_longer(cols = dplyr::everything(), 
            names_to = "condition", values_to = "value")

    tibble::tibble(score_values, set = contamination)
}


#' @name plotSampleScores
#' 
#' @title Plot the contamination scores of samples
#' 
#' @description
#' The function \code{plotSampleScores} creates a barplot or boxplot that 
#' visualizes the contamination score values obtained by 
#' \code{scoreSamples}.
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
#' @importFrom ggplot2 element_text ylab xlab
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
    if (!"value" %in% colnames(scores)) 
        stop("column 'value' not in 'scores'")
    if (!"condition" %in% colnames(scores)) 
        stop("column 'condition' not in 'scores'")
    if (!"set" %in% colnames(scores)) 
        stop("column 'set' not in 'scores'")
    
    if (any(duplicated(scores$condition))) {
        .method <- "permuteProteins"
    } else {
        .method <- "decoupleProteins"
    }
    
    ## factorize the condition column (used for x-axis)
    scores$condition <- factor(scores$condition, 
        levels = unique(scores$condition))
    
    ## do the actual plotting
    g <- ggplot2::ggplot(scores, 
        ggplot2::aes_string(y = "value", x = "condition"))
    
    ## depending on the function that was used to calculate the scores, plot 
    ## either a barplot (function decoupleProteins)or a boxplot 
    ## (function permuteProteins)
    if (.method == "decoupleProteins")
        g <- g + ggplot2::geom_bar(
            ggplot2::aes_string(fill = "set", col = "set"), 
            position = "dodge", stat = "identity", color="black", ) +
            ggplot2::ylab("s.d. to empirical null distribution")

    if (.method == "permuteProteins")
        g <- g + ggplot2::geom_boxplot(ggplot2::aes_string(x = "condition")) + 
            ggplot2::ylab("difference to markers")

    ## continue with changing the theme and return the plot
    g + ggplot2::theme_minimal() +
        ggplot2::xlab("samples") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, 
            angle = 90, hjust = 1))
}

