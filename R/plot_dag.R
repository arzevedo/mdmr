#' Plot graph or heatmap of Adjacency Matrix
#'
#' Plots the structure of a dynamic Bayesian network estimated by `mdm()` using either
#' a DAG (directed acyclic graph) representation or a heatmap of the adjacency matrix.
#'
#' @param mdm_object An object of class \code{mdm}, typically returned by the \code{mdm()} function.
#' @param node_labels Optional vector of node labels. If \code{NULL}, uses the column names of \code{mdm_object$data}, or defaults to \code{V1, V2, ..., Vn}.
#' @param type Type of plot to display. Either \code{"graph"} (default) for a directed graph with arrows or \code{"heatmap"} for a matrix view.
#' @param show_legend Logical. If \code{TRUE}, displays a legend for node/edge elements. Default is \code{FALSE}.
#' @param edge_color Color of the graph edges (arrows). Default is \code{"black"}.
#' @param node_color Fill color of the graph nodes. Default is \code{"steelblue"}.
#' @param label_color Color of the text labels inside the nodes. Default is \code{"white"}.
#' @param arrow_size Size of the graph arrow heads (in mm units). Default is \code{4}.
#'
#' @return A \code{ggplot} object representing either the graph or the heatmap from the adjacency matrix.
#'
#' @details
#' This function provides two modes of visualization:
#' \itemize{
#'   \item A \strong{heatmap} showing parent-child relations based on the adjacency matrix.
#'   \item A \strong{graph} using the \code{ggraph} and \code{igraph} packages.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_segment geom_text geom_point geom_label scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal coord_fixed theme element_blank scale_x_discrete theme_void
#' @importFrom grid arrow unit
#' @importFrom igraph graph_from_adjacency_matrix V
#' @importFrom ggraph ggraph geom_edge_link geom_node_label circle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- mdm(y)
#' plot_dag(res, type = "graph")
#' plot_dag(res, type = "heatmap")
#' }
plot_dag <- function(mdm_object, node_labels = NULL,
                      type = c("heatmap", "graph"),
                      show_legend = FALSE, edge_color = "black",
                      node_color = "steelblue", label_color = "white",
                      arrow_size = 4) {
  
  type <- match.arg(type)
  if (!inherits(mdm_object, "mdm")) stop("Input must be an object of class 'mdm'")
  adj_mat <- mdm_object$adj_mat
  n <- ncol(adj_mat)
  
  if (is.null(node_labels)) {
    node_labels <- if (!is.null(colnames(mdm_object$data))) colnames(mdm_object$data)
    else paste0("V", seq_len(n))
  }
  
  if (type == "heatmap") {
    df <- expand.grid(Child = node_labels, Parent = rev(node_labels),
                      stringsAsFactors = FALSE)
    df$value <- as.vector(t(adj_mat[rev(seq_len(n)), ]))
    df$Child <- factor(df$Child, levels = node_labels)
    df$Parent <- factor(df$Parent, levels = rev(node_labels))
    
    ggplot2::ggplot(df, ggplot2::aes(x = Child, y = Parent, fill = factor(value))) +
      ggplot2::geom_tile(color = "gray") +
      ggplot2::scale_fill_manual(values = setNames(c("white","red","grey80"),
                                                   c("0","1","NA")),
                                 name = NULL, breaks = c("1"), labels = c("Aresta")) +
      ggplot2::labs(x = "Filho", y = "Pai") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::theme(legend.position = "none",
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank()) +
      ggplot2::coord_fixed(expand = FALSE)
    
  } else {
    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed")
    igraph::V(g)$name <- node_labels
    
    ggraph::ggraph(g, layout = "sugiyama") +
      ggraph::geom_edge_link(
        arrow = grid::arrow(length = grid::unit(arrow_size, "mm"), type = "closed"),
        end_cap = ggraph::circle(3, "mm"),
        edge_width = 0.8,
        edge_colour = edge_color
      ) +
      ggraph::geom_node_label(
        ggplot2::aes(label = name),
        label.size = NA,
        fill = node_color,
        color = label_color,
        fontface = "bold",
        size = 4
      ) +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  }
}
