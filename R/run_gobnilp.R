#' Run GOBNILP Optimizer and Extract Adjacency Matrix
#'
#' Calls the GOBNILP binary on a precomputed `.scores` file and extracts the resulting 
#' adjacency matrix from the output.
#'
#' @param scores_path Character string. Path to the `.scores` file (Jaakkola format).
#' @param gobnilp_path Character string. Full path to the GOBNILP binary executable.
#' @param output_path Optional. Path to the output file for the adjacency matrix. If not provided, a temporary file will be used.
#'
#' @return A binary adjacency matrix representing the learned DAG.
#' @export
run_gobnilp <- function(scores_path, gobnilp_path,
                        output_path = tempfile(fileext = ".mat")) {
  
  setfile <- tempfile(fileext = ".set")
  writeLines(sprintf('gobnilp/outputfile/adjacencymatrix = "%s"', output_path), setfile)
  
  cmd <- sprintf('"%s" -f=jkl -v=2 -g="%s" "%s"', gobnilp_path, setfile, scores_path)
  system(cmd#, ignore.stdout = TRUE, ignore.stderr = TRUE
  )
  
  if (!file.exists(output_path)) stop("Adjacency matrix not generated.")
  mat <- as.matrix(read.table(output_path, header = FALSE))
  
  unlink(c(setfile, output_path))
  return(mat)
}