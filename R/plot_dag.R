#' Plot Static Structure of the Bayesian Network as a Heatmap
#'
#' Visualizes the adjacency matrix of the estimated Bayesian network from an MDM object.
#' Each tile indicates the presence (or absence) of a directed edge from a parent to a child node.
#'
#' @param mdm_object An object of class \code{"mdm"} returned by \code{\link{mdm}}.
#' @param node_labels Optional character vector to label nodes in the heatmap.
#'   If \code{NULL}, nodes are labeled as \code{"V1"}, \code{"V2"}, ..., \code{"VN"}.
#' @param fill_colors A character vector of three colors corresponding to the fill values for \code{0}, \code{1}, and \code{NA}.
#'
#' @return A \code{ggplot} object showing a heatmap of the adjacency matrix.
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme_minimal
#'   scale_x_discrete theme coord_fixed element_blank
#'
#' @examples
#' \dontrun{
#' # Estimate structure and plot adjacency heatmap
#' result <- mdm(data_input = y)
#' plot_dag(result)
#'
#' # Custom node labels
#' plot_dag(result, node_labels = c("X1", "X2", "X3"))
#' }
#' @export
plot_dag <- function(mdm_object, node_labels = NULL,
                             fill_colors = c("white", "red", "grey80")) {
  if (!inherits(mdm_object, "mdm")) stop("Input must be an object of class 'mdm'")
  if (!is.matrix(mdm_object$adj_mat)) stop("Input must be a matrix.")
  
  adj_mat <- mdm_object$adj_mat
  n <- ncol(adj_mat)
  if (is.null(node_labels)) {
    node_labels <- paste0("V", seq_len(n))
  }
  
  df <- expand.grid(Child = node_labels, Parent = rev(node_labels),
                    stringsAsFactors = FALSE)
  df$value <- as.vector(t(adj_mat[rev(seq_len(n)), ]))
  
  df$Child <- factor(df$Child, levels = node_labels)
  df$Parent <- factor(df$Parent, levels = rev(node_labels))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Child, y = Parent, fill = factor(value))) +
    ggplot2::geom_tile(color = "gray") +
    ggplot2::scale_fill_manual(values = setNames(fill_colors, c("0", "1", "NA")),
                               name = NULL,
                               breaks = c("1"),
                               labels = c("Edge")) +
    ggplot2::labs(x = "Child", y = "Parent") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(
      legend.position = "none",
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed(expand = FALSE)
  
  return(p)
}
