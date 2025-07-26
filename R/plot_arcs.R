#' Plot Dynamic Parameters of a Bayesian Network Over Time
#'
#' This function plots time series of either intercepts, connection strengths,
#' or both from an object of class \code{"mdm"} using posterior means and
#' credible/confidence intervals. The user can select between filtered or smoothed estimates.
#'
#' @param mdm_object An object of class \code{"mdm"} returned by \code{\link{mdm}}.
#' @param type Character string indicating which parameters to plot.
#'   Choices are:
#'   \describe{
#'     \item{"connections"}{Plots only dynamic edge parameters (e.g., V3->V1).}
#'     \item{"intercepts"}{Plots only intercepts (e.g., beta0 terms).}
#'     \item{"all"}{Plots both intercepts and connections.}
#'   }
#' @param distribution Whether to use filtered ("filt") or smoothed ("smoo") distributions.
#'   Default is "filt".
#' @param ci_level Numeric value between 0 and 1. Width of the confidence/credible interval.
#'   Default is 0.95.
#'
#' @return A \code{ggplot} object showing the posterior trajectories with uncertainty bounds.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon facet_wrap labs theme_minimal theme element_text
#' @importFrom stats qt
#'
#' @examples
#' \dontrun{
#' # Fit model
#' result <- mdm(data_input = y)
#'
#' # Plot dynamic connectivity (filtered)
#' plot_arcs(result, type = "connections")
#'
#' # Plot smoothed intercepts
#' plot_arcs(result, type = "intercepts", distribution = "smoo")
#' }
#' @export
plot_arcs <- function(mdm_object, type = c("connections", "intercepts", "all"), 
                                    distribution = "filt", ci_level = 0.95){
  if (!inherits(mdm_object, "mdm")) stop("Input must be of class 'mdm'")
  type <- match.arg(type)
  distribution <- match.arg(distribution, choices = c("filt", "smoo"))
  
  if(distribution == "filt"){
    mt_list <- mdm_object$Filt$mt
    Ct_list <- mdm_object$Filt$Ct
    nt_list <- mdm_object$Filt$nt
  } else {
    mt_list <- mdm_object$Smoo$smt
    Ct_list <- mdm_object$Smoo$sCt
    SE_list <- mdm_object$Smoo$SE
  }
  
  all_data <- list()
  
  for (node in seq_along(mt_list)){
    mt_node <- mt_list[[node]]
    Ct_node <- Ct_list[[node]]
    if(distribution == "filt") nt_node <- nt_list[[node]]
    if(distribution == "smoo") se_node <- SE_list[[node]]
    
    if(is.null(dim(mt_node))){
      mt_node <- matrix(mt_node, nrow = 1)
      rownames(mt_node) <- if(!is.null(names(mt_node))) names(mt_node) else paste0("beta0_", node)
    }
    if(is.null(rownames(mt_node))){
      rownames(mt_node) <- paste0("beta0_", node)
    }
    if(is.null(dim(Ct_node)) & distribution == "filt"){
      Ct_node <- array(Ct_node, dim = c(1, 1, length(mt_node)))
    }
    
    param_names <- rownames(mt_node)
    time <- seq_len(ncol(mt_node))
    
    for (param in seq_len(nrow(mt_node))) {
      name <- param_names[param]
      is_intercept  <- startsWith(name, "beta0")
      is_connection <- grepl("->", name)
      
      if ((type == "connections" && !is_connection) ||
          (type == "intercepts" && !is_intercept) ||
          (type == "all" && !(is_connection || is_intercept))) {
        next
      }
      
      mean_vals <- mt_node[param, ]
      var_vals  <- Ct_node[param, param, ]
      if (distribution == "filt"){
        tt <- qt((1 + ci_level) / 2, df = nt_node) * sqrt(var_vals)
      } else {
        if(length(dim(se_node))>2){
          tt <- se_node[,,]
        } else tt <- se_node[,param]
      }
      
      df <- data.frame(
        node      = paste0("Node ", node),
        parameter = name,
        time      = time,
        mean      = mean_vals,
        lower     = mean_vals - tt,
        upper     = mean_vals + tt,
        stringsAsFactors = FALSE
      )
      all_data[[length(all_data) + 1]] <- df
    }
  }
  
  if (length(all_data) == 0) stop("No parameters matched the selected type.")
  
  df_all <- do.call(rbind, all_data)
  df_all$facet_label <- paste0(df_all$node, ": ", df_all$parameter)
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = time, y = mean)) +
    ggplot2::geom_line(color = "navy") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "navy") +
    ggplot2::facet_wrap(~ facet_label, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("Posterior Estimates and %.0f%% CI (%s)", ci_level * 100, toupper(distribution)),
      x = "Time", y = "Posterior Estimate"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  
  return(p)
}
