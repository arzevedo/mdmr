% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdm_structure.R
\name{mdm_structure}
\alias{mdm_structure}
\title{Learn MDM Structure and Return Adjacency Matrix}
\usage{
mdm_structure(
  data_input,
  method = "hc",
  gobnilp_path = NULL,
  nbf = 15,
  delta = seq(0.5, 1, 0.01),
  subjects_length = 1,
  verbose = TRUE
)
}
\arguments{
\item{data_input}{Data frame or array. Input multivariate time series (possibly with subjects).}

\item{method}{Character. "hc" to use hill-climbing from `bnlearn`, or "ipa" to use the GOBNILP optimizer.}

\item{gobnilp_path}{Character. Path to the compiled GOBNILP binary (required if method = "gobnilp").}

\item{nbf}{Numeric. Starting time point for log-predictive likelihood.}

\item{delta}{Numeric vector. Sequence of discount factors.}

\item{subjects_length}{Numeric. Number of subjects in data.}

\item{verbose}{Logical. Whether to print progress.}
}
\value{
A binary adjacency matrix.
}
\description{
This function computes MDM scores and learns the structure of a dynamic Bayesian network,
returning its adjacency matrix using either `bnlearn::hc` or GOBNILP optimizer.
}
