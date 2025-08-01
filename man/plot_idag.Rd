% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_idag.R
\name{plot_idag}
\alias{plot_idag}
\title{Animate Dynamic Connectivity Heatmap Over Time}
\usage{
plot_idag(
  mdm_object,
  output_gif = "mdm.gif",
  fps = 10,
  width = 6,
  height = 6,
  dpi = 150
)
}
\arguments{
\item{mdm_object}{An object of class \code{"mdm"} as returned by \code{\link{mdm}}.}

\item{output_gif}{Name of the output file. Must end with ".gif". Default is \code{"mdm.gif"}.}

\item{fps}{Frames per second for the animation. Default is 10.}

\item{width}{Width (in inches) of each frame. Default is 6.}

\item{height}{Height (in inches) of each frame. Default is 6.}

\item{dpi}{Resolution (dots per inch) for saved frames. Default is 150.}
}
\value{
A \code{magick-image} object containing the animated heatmap.
}
\description{
Creates and saves a GIF showing the evolution of dynamic connectivity
in a Bayesian network learned with the \code{\link{mdm}} function.
}
\details{
The animation reflects time-varying intensity of edges in the network,
based on the dynamic parameters.
}
\examples{
\dontrun{
# Estimate MDM model
result <- mdm(data_input = y)

# Create and save the animated heatmap
plot_idag(result, fps = 12, output_gif = "network_evolution.gif")
}
}
