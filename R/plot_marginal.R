#' Plot Dynamic Influence of Parents on a Target Node
#'
#' Visualizes the dynamic influence of each parent node on the target node
#' over time by computing the time-varying contribution of each parent series
#' weighted by its dynamic coefficient.
#'
#' @param mdm_object An object of class \code{"mdm"} as returned by \code{\link{mdm}}.
#' @param target_node Integer index specifying the target node (child).
#' @param distribution Whether to use "filt" (filtered) or "smoo" (smoothed) posterior estimates.
#'   Default is \code{"filt"}.
#' @param scale_series Logical. If \code{TRUE}, all time series (observed and parental contributions)
#'   are standardized (mean zero, unit variance). Default is \code{FALSE}.
#'
#' @details
#' The contribution of a parent to the target node at time \eqn{t} is calculated as:
#' \deqn{
#' \text{Contribution}_{i \rightarrow j}(t) = \beta_{i \rightarrow j}(t) \cdot Y_i(t)
#' }
#' where \eqn{\beta_{i \rightarrow j}(t)} is the time-varying coefficient estimated from the dynamic model,
#' and \eqn{Y_i(t)} is the observed time series of parent node \eqn{i}.
#'
#' These contributions are plotted as colored lines alongside the observed series for the target node (in black),
#' enabling interpretation of how the influence from each parent evolves over time.
#'
#' @return A \code{ggplot} object showing the dynamic contributions and the observed series.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal scale_color_brewer
#'
#' @examples
#' \dontrun{
#' # Estimate model
#' result <- mdm(data_input = y)
#'
#' # Plot parental influence (no scaling)
#' plot_marginal(result, Y = as.matrix(y), target_node = 2, scale_series = FALSE)
#'
#' # Standardize series before plotting
#' plot_marginal(result, Y = as.matrix(y), target_node = 2, scale_series = TRUE)
#' }
#' @export
plot_marginal <- function(mdm_object, target_node, distribution = "filt",
                                   scale_series = FALSE) {
  if (!inherits(mdm_object, "mdm")) stop("mdm_object must be of class 'mdm'")
  if(is.null(target_node)){
    stop("Select a node to verify the marginal contribution")
    }
  if(is.null(mdm_object$data)) {
    stop("mdm_object must contain a 'data' matrix as 'data' field.")
  }
  Y <- mdm_object$data
  distribution <- match.arg(distribution, choices = c("filt", "smoo"))
  
  if (distribution == "filt") {
    beta_list <- mdm_object$Filt$mt
  } else {
    beta_list <- mdm_object$Smoo$smt
  }
  
  beta_node <- beta_list[[target_node]]
  if (is.null(dim(beta_node))) stop("No dynamic parameters found for this node.")
  
  param_names <- rownames(beta_node)
  Time <- ncol(beta_node)
  time <- seq_len(Time)
  
  # Identify dynamic connections (parent → target_node)
  parent_rows <- which(grepl("->", param_names))
  if (length(parent_rows) == 0) stop("No parent connections found for this node.")
  
  child_series <- Y[, target_node]
  if (scale_series) child_series <- scale(child_series)[,1]
  
  df_list <- list()
  for (row in parent_rows){
    name <- param_names[row]
    parent <- strsplit(name, "->")[[1]][1]
    parent_idx <- which(colnames(Y) == parent)
    
    if (length(parent_idx) != 1){
      warning(sprintf("Parent %s not found in column names of Y — skipping", parent))
      next
    }
    
    parent_series <- Y[, parent_idx]
    if (scale_series) parent_series <- scale(parent_series)[,1]
    
    beta_vals <- beta_node[row, ]
    contribution <- beta_vals * parent_series
    
    df_list[[length(df_list) + 1]] <- data.frame(
      time = time,
      contribution = contribution,
      parent = parent
    )
  }
  
  df_contrib <- do.call(rbind, df_list)
  df_child <- data.frame(
    time = time,
    value = child_series,
    series = "Observed"
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_contrib, ggplot2::aes(x = time, y = contribution, color = parent),
                       alpha = 0.8, linewidth=.9) +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::geom_line(data = df_child,
                       ggplot2::aes(x = time, y = value), color = "black", linewidth = .5, linetype=1) +
    ggplot2::labs(
      title = sprintf("Dynamic Parental Influence on Node %d (%s)", target_node, toupper(distribution)),
      y = "Contribution / Observed",
      x = "Time", color = "Parent"
    ) +
    ggplot2::theme_minimal(base_size = 13)
  
  return(p)
}
