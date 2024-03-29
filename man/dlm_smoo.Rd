% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlm_smoo.R
\name{dlm_smoo}
\alias{dlm_smoo}
\title{Smooth dynamic linear models (DLMs)}
\usage{
dlm_smoo(mt, Ct, Rt, nt, dt, Gt = 0)
}
\arguments{
\item{mt}{The matrix of posterior mean with dimension p X T}

\item{Ct}{The squared matrix of posterior variance with dimension p X p X T}

\item{Rt}{The squared matrix of prior variance with dimension p X p X T}

\item{nt}{The vector of prior hypermarameters of precision phi with length T}

\item{dt}{The vector of prior hypermarameters of precision phi with length T}

\item{Gt}{he matrix of state equation with dimension: p X p X T. The default is identity matrix block}
}
\description{
DLM with Smoothing for unknown observational and state variances
}
