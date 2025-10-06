#' Learn MDM Structure and Return Adjacency Matrix
#'
#' @description
#' This function computes MDM scores and learns the structure of a dynamic Bayesian network,
#' returning its adjacency matrix using either `bnlearn::`(hc/tabu/mmhc/h2pc/rsmax2) with a custom 
#' decomposable or GOBNILP ("ipa") optimizer.
#'
#' @param data_input Data frame. Input multivariate time series (possibly with subjects).
#' @param method A \code{character} string indicating the method to use for structure learning. Options are:
#'   \describe{
#'     \item{"hc"}{Hill-climbing structure learning using the `bnlearn` package.}
#'     \item{"tabu"}{Tabu search greedy search learning using the `bnlearn` package.}
#'     \item{"h2pc"}{Hybrid HPC (H2PC) learning using the `bnlearn` package.}
#'     \item{"mmhc"}{Max-Min Hill Climbing (MMHC) learning using the `bnlearn` package.}
#'     \item{"rsmax2"}{2-phase Restricted Maximization (RSMAX2) learning using the `bnlearn` package.}
#'     \item{"ipa"}{Integer Programming Approach using GOBNILP with Jaakkola scoring.}
#'   }
#' @param gobnilp_path Character. Path to the compiled GOBNILP binary (required if method = "gobnilp").
#' @param nbf Numeric. Starting time point for log-predictive likelihood.
#' @param delta Numeric vector. Sequence of discount factors.
#' @param subjects_length Numeric. Number of subjects in data.
#' @param verbose Logical. Whether to print progress.
#' @return A binary adjacency matrix.
mdm_structure <- function(data_input,
                          method = "hc",
                          gobnilp_path = NULL,
                          nbf = 15,
                          delta = seq(0.5, 1.0, 0.01),
                          subjects_length = 1,
                          verbose = TRUE) {
  
  method <- match.arg(method, choices = c(
    "hc", "tabu", "mmhc", "h2pc", "rsmax2",
    "ipa")
  )
  
  # bnlearn >= 5.0 uses "custom-score"; older versions used "custom".
  bn_ver <- as.character(utils::packageVersion("bnlearn"))
  if(package_version(packageVersion("bnlearn"))>"4.9"){
    arg_score_name <- "custom-score"
  }else{
    arg_score_name <- "custom"
  }
  
  if (method %in% c("hc", "tabu", "mmhc", "h2pc", "rsmax2")) {
    if (verbose) message(sprintf("Running bnlearn::%s() with MDM custom score...", method))
    
    if (is.null(colnames(data_input))) {
      colnames(data_input) <- paste0("V", seq_len(ncol(data_input)))
    }
    dat <- data.frame(data_input)
    
    bn_obj <- switch(
      method,
      "hc" = bnlearn::hc(
        x     = dat,
        score = arg_score_name,
        fun   = mdm_score_bn,
        args  = list(nbf = nbf, method = "Brent", call = FALSE)
      ),
      "tabu" = bnlearn::tabu(
        x     = dat,
        score = arg_score_name,
        fun   = mdm_score_bn,
        args  = list(nbf = nbf, method = "Brent", call = FALSE)
      ),
      "mmhc" = bnlearn::mmhc(
        x = dat,
        restrict.args = list(), 
        maximize.args = list(
          score = arg_score_name,
          fun   = mdm_score_bn,
          args  = list(nbf = nbf, method = "Brent", call = FALSE)
        )
      ),
      "h2pc" = bnlearn::h2pc(
        x = dat,
        restrict.args = list(),
        maximize.args = list(
          score = arg_score_name,
          fun   = mdm_score_bn,
          args  = list(nbf = nbf, method = "Brent", call = FALSE)
        )
      ),
      "rsmax2" = bnlearn::rsmax2(
        x = dat,
        restrict.args = list(),
        maximize.args = list(
          score = arg_score_name,
          fun   = mdm_score_bn,
          args  = list(nbf = nbf, method = "Brent", call = FALSE)
        )
      )
    )
    
    N <- ncol(dat)
    var_names <- colnames(dat)
    adj <- matrix(0L, nrow = N, ncol = N, dimnames = list(var_names, var_names))
    
    if (!is.null(bn_obj$arcs) && nrow(bn_obj$arcs) > 0L) {
      for (z in seq_len(nrow(bn_obj$arcs))) {
        from <- bn_obj$arcs[z, 1]
        to   <- bn_obj$arcs[z, 2]
        adj[from, to] <- 1L
      }
    }
    
  } else if (method == "ipa") {
    if (is.null(gobnilp_path)) stop("You must provide gobnilp_path when using method = 'ipa'.")
    if (verbose) message("Computing MDM scores for GOBNILP (IPA) pipeline...")
    colnames(data_input) <- NULL
    
    score_result <- mdm_score(data_input, nbf = nbf, delta = delta,
                              GOLB_print = TRUE, subjects_length = subjects_length)
    
    aux_score <- file.info(list.files(pattern = "^mdm_score_.*", full.names = T))
    score_file <- rownames(aux_score)[which.max(aux_score$mtime)]
    
    if (verbose) message("Running GOBNILP...")
    adj <- run_gobnilp(scores_path = score_file, gobnilp_path = gobnilp_path)
    rownames(adj) <- colnames(adj)
    unlink(score_file)
  }
  
  return(adj)
}
