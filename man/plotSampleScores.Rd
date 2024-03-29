% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantifyViaDecoupleR.R
\name{plotSampleScores}
\alias{plotSampleScores}
\title{Plot the contamination scores of samples}
\usage{
plotSampleScores(scores)
}
\arguments{
\item{scores}{tibble}
}
\value{
\code{gg}
}
\description{
The function \code{plotSampleScores} creates a barplot or boxplot that 
visualizes the contamination score values obtained by 
\code{scoreSamples}.
}
\details{
The function takes the argument \code{scores} that can be created via 
\code{scoreSamples}.
}
\examples{
library(SummarizedExperiment)
library(tibble)

f <- system.file("protein_datasets/tanzer2020.RDS", 
    package = "apoptosisQuantification")
tanzer2020 <- readRDS(f)
prot <- assay(tanzer2020)

## use the intensities to calculate scores
scores <- scoreSamples(values = prot, contamination = "apoptosis",fc = 2, n = 1, n_perm = 100)
scores$treatment <- tanzer2020$treatment
plotSampleScores(scores = scores) +
    ggplot2::facet_wrap(~ treatment, scales = "free_x") +
    ggplot2::theme(legend.position = "none")

## use the intensities to calculate scores
scores <- scoreSamples(values = prot, contamination = "necroptosis", n_perm = 100)
scores$treatment <- tanzer2020$treatment
plotSampleScores(scores = scores) +
    ggplot2::facet_wrap(~ treatment, scales = "free_x") +
    ggplot2::theme(legend.position = "none")
}
