library(SummarizedExperiment)
## tanzer2020_secretome
## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
## download: u937_apoptosis_9h12-5
tanzer2020 <- MatrixQCvis::maxQuant("u937_apoptosis_9h12-5/proteinGroups.txt", 
    type = "txt", intensity = "LFQ")

## data transformation
tanzer2020 <- tanzer2020[-grep(rownames(tanzer2020), pattern = "REV|CON__"), ]
a_t <- MatrixQCvis::transformAssay(a = assay(tanzer2020), method = "vsn")
a_i <- MatrixQCvis::imputeAssay(a = a_t, method = "MinDet")
tanzer2020 <- MatrixQCvis:::updateSE(se = tanzer2020, assay = a_i)

## translate the proteins
library("org.Hs.eg.db")
uniprots <- Rkeys(org.Hs.egUNIPROT)
dict <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "SYMBOL", "UNIPROT")
uniprot <- dict$UNIPROT
symbol <- dict$SYMBOL
names(symbol) <- uniprot

## rename proteins
rownames(tanzer2020) <- unlist(lapply(rownames(tanzer2020), 
    function(x) paste(symbol[strsplit(x, split = ";")[[1]]], collapse = ";")))

## remove the features that have no corresponding Symbol
exclude_pattern <- c("NA", "NA;NA", "NA;NA;NA", "NA;NA;NA;NA")
tanzer2020 <- tanzer2020[!rownames(tanzer2020) %in% exclude_pattern, ]

## only take the control and TNF_SM (apoptosis-induced) samples
samps <- grep(colnames(tanzer2020), pattern = "TNF_SM_|DMSO")
tanzer2020 <- tanzer2020[, samps]

## remove un-paired samples
samps <- grep(colnames(tanzer2020), pattern = "DMSO_9h|DMSO_12[.]5h")
tanzer2020 <- tanzer2020[, -samps]

## rename the columns
colnames(tanzer2020) <- unlist(lapply(
    strsplit(colnames(tanzer2020), split = "LFQ[.]intensity[.]"), "[", 2))

## modify colData
tanzer2020$name <- tanzer2020$name_cut <- colnames(tanzer2020)
tanzer2020$treatment <- ifelse(grepl(tanzer2020$name, pattern = "DMSO"), 
    "DMSO", "TNF_SM")
tanzer2020$time <- unlist(lapply(strsplit(tanzer2020$name, split = "_"), 
    function(name_i) name_i[length(name_i) - 1]))

## save the file
saveRDS(tanzer2020, file = "tanzer2020_secretome.RDS")

## tanzer2020
## source: https://www.ebi.ac.uk/pride/archive/projects/PXD014966
## download: u937_necro_a_cell
tanzer2020 <- MatrixQCvis::maxQuant("u937_necro_a_cell/proteinGroups.txt", 
    type = "txt", intensity = "LFQ")

## data transformation
tanzer2020 <- tanzer2020[-grep(rownames(tanzer2020), pattern = "REV|CON__"), ]
a_t <- MatrixQCvis::transformAssay(a = assay(tanzer2020), method = "vsn")
a_i <- MatrixQCvis::imputeAssay(a = a_t, method = "MinDet")
tanzer2020 <- MatrixQCvis:::updateSE(se = tanzer2020, assay = a_i)

## translate the proteins
library("org.Hs.eg.db")
uniprots <- Rkeys(org.Hs.egUNIPROT)
dict <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "SYMBOL", "UNIPROT")
uniprot <- dict$UNIPROT
symbol <- dict$SYMBOL
names(symbol) <- uniprot

## rename proteins
rownames(tanzer2020) <- unlist(lapply(rownames(tanzer2020), 
    function(x) paste(symbol[strsplit(x, split = ";")[[1]]], collapse = ";")))

## remove the features that have no corresponding Symbol
exclude_pattern <- c("NA", "NA;NA", "NA;NA;NA", "NA;NA;NA;NA")
tanzer2020 <- tanzer2020[!rownames(tanzer2020) %in% exclude_pattern, ]

## only take the control and TNF_SM (apoptosis-induced) samples
#samps <- grep(colnames(tanzer2020), pattern = "TNF_SM_|DMSO")
#tanzer2020 <- tanzer2020[, samps]

## rename the columns
colnames(tanzer2020) <- unlist(lapply(
    strsplit(colnames(tanzer2020), split = "LFQ[.]intensity[.]"), "[", 2))

## modify colData
tanzer2020$name <- colnames(tanzer2020)
tanzer2020$treatment <- unlist(lapply(
    strsplit(tanzer2020$name, split = "_1h|_3h|_5h|_7h"), "[", 1))
tanzer2020$time <- unlist(lapply(strsplit(tanzer2020$name, split = "_"), 
                                 function(name_i) name_i[length(name_i) - 1]))

## save the file
saveRDS(tanzer2020, file = "tanzer2020.RDS")
c