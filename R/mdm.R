#' Estimate Bayesian Network Structure and Dynamic Parameters (MDM)
#'
#' This is the main function of the MDMr package. It estimates the structure of a Bayesian Network from multivariate time series data
#' and calculates dynamic parameters using both filtering and smoothing techniques. The structure can be inferred using either the 
#' hill-climbing algorithm from the `bnlearn` package or using `IPA` via a GGOBNILP and Jaakkola local score file interface.
#'
#' @param data_input A \code{data.frame} of observed data. Rows represent time points, and columns represent nodes. Must be complete (no missing values).
#' @param method A \code{character} string indicating the method to use for structure learning. Options are:
#'   \describe{
#'     \item{"hc"}{Hill-climbing structure learning using the `bnlearn` package.}
#'     \item{"ipa"}{Integer Programming Approach using GOBNILP with Jaakkola scoring.}
#'   }
#' @param ... Additional arguments passed to \code{\link{mdm_structure}}, e.g., \code{gobnilp_path}.
#'
#' @return An object of class \code{"mdm"}, which is a list with the following components:
#' \describe{
#'   \item{adj_mat}{Adjacency matrix representing the learned DAG structure.}
#'   \item{DF}{A list with discount factor estimation results.}
#'   \item{Filt}{A list of filtered dynamic parameters.}
#'   \item{Smoo}{A list of smoothed dynamic parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run with hill-climbing
#' result_hc <- mdm(data_input = y)
#'
#' # Run with GOBNILP
#' result_ipa <- mdm(data_input = y, method = "ipa",
#'                   gobnilp_path = "/path/to/gobnilp")
#' }
#' @export
mdm <- function(data_input, method = "hc", ...) {
  if (!is.data.frame(data_input)) {
    stop("data_input must be a dataframe of observed data [T × N]")
  }
  method <- match.arg(method, choices = c("hc", "ipa"))
  
  # Structure (DAG)
  adj_mat <- mdm_structure(data_input = data_input, method = method, ...)
  # Discount Factor
  DF1 <- CDELT(as.matrix(data_input), adj_mat)
  # Filtered Distributions
  Filt <- mdm_filt(as.matrix(data_input), adj_mat, DF1$DF_hat)
  # Smoothed Distributions
  Smoo <- mdm_smoo(Filt$mt, Filt$Ct, Filt$Rt, Filt$nt, Filt$dt)
  
  # Output
  output <- list(
    adj_mat = adj_mat,
    DF = DF1,
    Filt = Filt,
    Smoo = Smoo,
    data = as.matrix(data_input)
  )
  class(output) <- "mdm"
  return(output)
}
