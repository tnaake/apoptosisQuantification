## load libraries
library("SummarizedExperiment")
library("org.Hs.eg.db")

## define helper functions

#' @name renameColumns
#' 
#' @description 
#' Remove \code{"LFQ.intensity."} from \code{colnames(se)}.
#' 
#' @param se SummarizedExperiment
#' 
#' @return character
renameColumns <- function(se) {
    unlist(lapply(strsplit(colnames(se), split = "LFQ[.]intensity[.]"), "[", 2))
}

#' @name modifyColData
#' 
#' @description 
#' Update the column \code{name} of \code{colData(se)} by \code{colnames(se)}.
#' 
#' Create the column \code{treatment} in \code{colData(se)} by taking the
#' respective treatment from \code{colData(se)$name}.
#' 
#' Create the column \code{time} in \code{colData(se)} by taking the 
#' respective time from \code{colData(se)$name}.
#' 
#' Create the column \code{treatment_time} in \code{colData(se)} by combining
#' the columns \code{treatment} and \code{time} of \code{colData(se)}.
#' 
#' @param se SummarizedExperiment
#' 
#' @return DataFrame
modifyColData <- function(se){
    cD <- colData(se)
    cD$name <- colnames(se)
    cD$treatment <- unlist(lapply(strsplit(cD$name, split = "_1h|_3h|_5h|_7h"), 
        "[", 1))
    cD$time <- unlist(lapply(strsplit(cD$name, split = "_"), 
        function(name_i) name_i[length(name_i) - 1]))
    cD$treatment_time <- paste(cD$treatment, cD$time, sep = "_")
    cD
}

#' @name modifyAssay
#' 
#' @description 
#' Filter the data set for at least 75% of valid values among three or 
#' four biological replicates in at least one condition.
#' 
#' Remove the proteins that are REV or COn sequences.
#' 
#' log-transform the data and impute by MinDet.
#' 
#' @param se SummarizedExperiment
#' 
#' @return SummarizedExperiment
modifyAssay <- function(se) {
    ## Quantified proteins are filtered for at least 75% of valid values among 
    ## three or four biological replicates in at least one condition. 
    measuredFeatures <- MatrixQCvis::measuredCategory(se = se, 
        measured = TRUE, category = "treatment_time")
    stopifnot(all(
        colnames(measuredFeatures[, -1]) == names(table(se$treatment_time))))
    measuredFeatures75 <- apply(measuredFeatures[, -1], 1, 
        function(x) any(as.vector(x) / table(se$treatment_time) >= 0.75))
    se <- se[measuredFeatures75, ]
    
    ## remove the proteins that are REV or CON sequences
    se <- se[-grep(rownames(se), pattern = "REV|CON__"), ]
    
    ## log-transform the data and impute by MinDet
    a_i <- MatrixQCvis::transformAssay(a = assay(se), method = "log") |>
        MatrixQCvis::imputeAssay(method = "MinDet")
    
    ## return the updates SummarizedExperiment object
    MatrixQCvis:::updateSE(se = se, assay = a_i)
}

#' @name simplifyRowData
#' 
#' @description 
#' Remove all columns in \code{rowData(se)} expect \code{feature}.
#'
#' @param se SummarizedExperiment
#'
#' @return SummarizedExperiment
simplifyRowData <- function(se) {
    ## remove unnecessary information in rowData
    rowData(se) <- rowData(se) |> 
        as.data.frame() |>
        dplyr::select(feature) |>
        S4Vectors::DataFrame()
    se
}

#' @name renameProteins
#' 
#' @description 
#' Translate the proteins from UNIPROT to SYMBOL ids.
#' 
#' Remove the proteins that have after translating the ids only NA values.
#' 
#' Remove the proteins with duplicated protein names.
#' 
#' @param se SummarizedExperiment
#' 
#' @return SummarizedExperiment
renameProteins <- function(se) {
    
    ## translate proteins
    uniprots <- Rkeys(org.Hs.egUNIPROT)
    dict <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "SYMBOL", "UNIPROT")
    uniprot <- dict$UNIPROT
    symbol <- dict$SYMBOL
    names(symbol) <- uniprot
    
    ## rename proteins
    featureNames <- rownames(se)
    featureNames <- lapply(featureNames, 
        function(name_i) strsplit(name_i, split = ";")[[1]])
    featureNames <- lapply(featureNames,
        function(name_i) paste(symbol[name_i], collapse = ";"))
    rownames(se) <- unlist(featureNames)
    
    ## remove the features that have no corresponding Symbol
    exclude_pattern <- c("NA", "NA;NA", "NA;NA;NA", "NA;NA;NA;NA", 
        "NA;NA;NA;NA;NA", "NA;NA;NA;NA;NA;NA", "NA;NA;NA;NA;NA;NA;NA",
        "NA;NA;NA;NA;NA;NA;NA;NA")
    se <- se[!rownames(se) %in% exclude_pattern, ]
    
    ## remove the duplicated proteins
    se[!duplicated(rownames(se)), ]
}

################################################################################
## tanzer2020_secretome
## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
## download: u937_apoptosis_9h12-5
tanzer2020 <- MatrixQCvis::maxQuant("u937_apoptosis_9h12-5/proteinGroups.txt", 
    type = "txt", intensity = "LFQ")

## modify colData
## rename the columns
colnames(tanzer2020) <- renameColumns(se = tanzer2020)
colData(tanzer2020) <- modifyColData(se = tanzer2020)

## data transformation
tanzer2020 <- modifyAssay(se = tanzer2020)

## only take the control and TNF_SM (apoptosis-induced) samples
samps <- grep(colnames(tanzer2020), pattern = "TNF_SM_|DMSO")
tanzer2020 <- tanzer2020[, samps]

## remove un-paired samples
samps <- grep(colnames(tanzer2020), pattern = "DMSO_9h|DMSO_12[.]5h")
tanzer2020 <- tanzer2020[, -samps]

## simplify rowData
tanzer2020 <- simplifyRowData(tanzer2020)

## rename proteins
tanzer2020 <- renameProteins(se = tanzer2020)

## save the file
saveRDS(tanzer2020, file = "tanzer2020_secretome.RDS")

## tanzer2020_necro_a_cell
## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
## download: u937_necro_a_cell
tanzer2020 <- MatrixQCvis::maxQuant("u937_necro_a_cell/proteinGroups.txt", 
    type = "txt", intensity = "LFQ")

## modify colData
## rename the columns
colnames(tanzer2020) <- renameColumns(se = tanzer2020)
colData(tanzer2020) <- modifyColData(se = tanzer2020)

## data transformation
tanzer2020 <- modifyAssay(se = tanzer2020)

## remove all columns except DMSO_, TNF_SM_, TNF_SM_Idun
tanzer2020 <- tanzer2020[, grepl(colnames(tanzer2020), pattern = "DMSO_[1|3|5|7]|TNF_SM_[1|3|5|7]|TNF_SM_Idun_[1|3|5|7]")]

## simplify rowData
tanzer2020 <- simplifyRowData(tanzer2020)

## rename proteins
tanzer2020 <- renameProteins(se = tanzer2020)

## save the file
saveRDS(tanzer2020, file = "tanzer2020_necro_a_cell.RDS")

################################################################################
## tanzer2020
## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
## download: u937_necro_apoptosis_super_search
tanzer2020 <- MatrixQCvis::maxQuant("u937_necro_apoptosis_super_search/proteinGroups.txt", 
        type = "txt", intensity = "LFQ")

## modify colData
## rename the columns
colnames(tanzer2020) <- renameColumns(se = tanzer2020)
colData(tanzer2020) <- modifyColData(se = tanzer2020)

## data transformation
tanzer2020 <- modifyAssay(se = tanzer2020)

## remove all columns except DMSO_, TNF_SM_, TNF_SM_Idun
tanzer2020 <- tanzer2020[, grepl(colnames(tanzer2020), pattern = "DMSO_[1|3|5|7]|TNF_SM_[1|3|5|7]|TNF_SM_Idun_[1|3|5|7]")]

## simplify rowData
tanzer2020 <- simplifyRowData(tanzer2020)

## rename proteins
tanzer2020 <- renameProteins(se = tanzer2020)

## save the file
saveRDS(tanzer2020, file = "tanzer2020.RDS")
