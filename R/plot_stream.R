#' Stream Plot of Parental Influence Over Time for a Target Node
#'
#' Visualizes the dynamic influence of each parent on a specific child node in a Bayesian dynamic network,
#' using streamgraph-style visualization based on posterior mean values.
#'
#' This function is useful for analyzing how the strength of influence from different parent nodes varies
#' across time for a specific target node.
#'
#' @param mdm_object An object of class \code{"mdm"} as returned by \code{\link{mdm}}.
#' @param child_node Integer index of the target node (i.e., the child whose parents' effects are shown).
#' @param distribution Whether to use filtered ("filt") or smoothed ("smoo") posterior estimates.
#'   Default is \code{"filt"}.
#'
#' @return A \code{ggplot} object containing the stream plot of parental influences.
#'
#' @importFrom ggplot2 ggplot aes labs scale_fill_brewer theme_minimal theme element_blank
#' @importFrom ggstream geom_stream
#'
#' @examples
#' \dontrun{
#' # Fit model
#' result <- mdm(data_input = y)
#'
#' # Plot stream graph of influence on node 2
#' plot_stream(result, child_node = 2)
#'
#' }
#' @export
plot_stream <- function(mdm_object, child_node, distribution = "filt"){
  if (!inherits(mdm_object, "mdm")) stop("Input must be of class 'mdm'")
  distribution <- match.arg(distribution, choices = c("filt", "smoo"))
  
  if (distribution == "filt") {
    mt_list <- mdm_object$Filt$mt
  } else {
    mt_list <- mdm_object$Smoo$smt
  }
  
  mt_node <- mt_list[[child_node]]
  
  param_names <- rownames(mt_node)
  
  if(!is.matrix(mt_node)){
    warning("Only dynamic paramater for this node is the Intercept. Defaulting to plot all intercepts")
    plot_arcs(mdm_object=mdm_object,type = "intercept")
  }else{
    time <- seq_len(ncol(mt_node))
    
    df_list <- list()
    
    for (param in seq_len(nrow(mt_node))) {
      name <- param_names[param]
      #if (!grepl("->", name)) next  # Skip intercepts
      
      parent <- strsplit(name, "->")[[1]][1]
      
      df <- data.frame(
        time   = time,
        parent = parent,
        value  = mt_node[param, ],
        stringsAsFactors = FALSE
      )
      df_list[[length(df_list) + 1]] <- df
    }
    
    if (length(df_list) == 0) stop("No parent connections found for this node.")
    
    df_all <- do.call(rbind, df_list)
    
    if (length(unique(df_all$parent)) == 1) {
      warning(sprintf(
        "Only one parent connection found for node %d - stream comparison may not be meaningful.",
        child_node))
    }
    
    p <- ggplot2::ggplot(df_all, ggplot2::aes(x = time, y = value, fill = parent)) +
      ggstream::geom_stream() +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(
        title = sprintf("Evolution of Parental Influence on Node %d (%s)", child_node, toupper(distribution)),
        x = "Time", y = "Posterior Mean", fill = "Parent"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
    
    return(p)
  }
}
