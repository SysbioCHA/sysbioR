% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LimmaDEG.R
\name{LimmaDEG}
\alias{LimmaDEG}
\title{DEG analaysis using Limma}
\usage{
LimmaDEG(motherMatrix, case, control)
}
\arguments{
\item{motherMatrix}{row name should be probe or gene and col name should be sample names.
mothermatrix colnames: treat_A, treat_B, treat_C, no_treat_A, no_treat_B, no_treat_C}

\item{case}{pattern of case samples}

\item{control}{pattern of control samples}
}
\value{
res DEG results.
}
\description{
This function performs DEG analysis using
}
