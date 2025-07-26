#' Animate Dynamic Connectivity Heatmap Over Time
#'
#' Creates and saves a GIF showing the evolution of dynamic connectivity
#' in a Bayesian network learned with the \code{\link{mdm}} function.
#'
#' The animation reflects time-varying intensity of edges in the network,
#' based on the dynamic parameters.
#'
#' @param mdm_object An object of class \code{"mdm"} as returned by \code{\link{mdm}}.
#' @param output_gif Name of the output file. Must end with ".gif". Default is \code{"mdm.gif"}.
#' @param fps Frames per second for the animation. Default is 10.
#' @param width Width (in inches) of each frame. Default is 6.
#' @param height Height (in inches) of each frame. Default is 6.
#' @param dpi Resolution (dots per inch) for saved frames. Default is 150.
#'
#' @return A \code{magick-image} object containing the animated heatmap.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal
#'   scale_x_discrete scale_y_discrete theme coord_fixed ggsave
#' @importFrom magick image_read image_animate image_write
#'
#' @examples
#' \dontrun{
#' # Estimate MDM model
#' result <- mdm(data_input = y)
#'
#' # Create and save the animated heatmap
#' plot_heatmap_animation(result, fps = 12, output_gif = "network_evolution.gif")
#' }
#' @export
plot_heatmap_animation <- function(mdm_object,
                                   output_gif = "mdm.gif",
                                   fps = 10,
                                   width = 6,
                                   height = 6,
                                   dpi = 150) {
  if (!inherits(mdm_object, "mdm")) stop("Input must be an object of class 'mdm'")
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("The 'magick' package is required. Install it with install.packages('magick')")
  }
  
  # Inlined data preparation
  tidy_all <- list()
  for (i in seq_along(mdm_object$Smoo$smt)) {
    mat <- mdm_object$Smoo$smt[[i]]
    if (is.null(rownames(mat))) rownames(mat) <- paste0("beta0_V", i)
    
    n_rows <- nrow(mat)
    n_cols <- ncol(mat)
    components <- rownames(mat)
    times <- seq_len(n_cols)
    
    df <- data.frame(
      node      = rep(i, each = n_rows * n_cols),
      component = rep(components, times = n_cols),
      time      = rep(times, each = n_rows),
      value     = as.vector(mat),
      stringsAsFactors = FALSE
    )
    
    df$parent <- sapply(df$component, function(comp) {
      if (grepl("->", comp)) strsplit(comp, "->")[[1]][1] else "beta0"
    })
    
    df$child <- sapply(df$component, function(comp) {
      if (grepl("->", comp)) {
        strsplit(comp, "->")[[1]][2]
      } else if (grepl("_", comp)) {
        sub(".*_", "", comp)
      } else {
        NA_character_
      }
    })
    
    tidy_all[[i]] <- df
  }
  
  df_long <- do.call(rbind, tidy_all)
  rownames(df_long) <- NULL
  
  # Remove intercepts
  df_long <- df_long[!grepl("beta0", df_long$parent), ]
  
  # Consistent axis levels
  all_labels <- sort(unique(c(df_long$parent, df_long$child)))
  df_long$parent <- factor(df_long$parent, levels = rev(all_labels))
  df_long$child  <- factor(df_long$child,  levels = all_labels)
  
  global_min <- min(df_long$value, na.rm = TRUE)
  global_max <- max(df_long$value, na.rm = TRUE)
  
  temp_dir <- tempfile(pattern = "frames_")
  dir.create(temp_dir)
  time_points <- sort(unique(df_long$time))
  
  for (i in seq_along(time_points)) {
    t <- time_points[i]
    plot_df <- df_long[df_long$time == t, ]
    plot_df$parent <- factor(plot_df$parent, levels = rev(all_labels))
    plot_df$child  <- factor(plot_df$child,  levels = all_labels)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = child, y = parent, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        name = "Connection Strength",
        limits = c(global_min, global_max)
      ) +
      ggplot2::scale_x_discrete(drop = FALSE, position = "top") +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::labs(x = "Child", y = "Parent", title = paste("Time:", t)) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::coord_fixed(expand = FALSE)
    
    frame_path <- file.path(temp_dir, sprintf("frame_%03d.png", i))
    ggplot2::ggsave(frame_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  imgs <- magick::image_read(list.files(temp_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE))
  anim <- magick::image_animate(imgs, fps = fps, optimize = TRUE)
  magick::image_write(anim, path = output_gif)
  gc()
  
  unlink(temp_dir, recursive = TRUE)
  return(anim)
}
