#' Learn MDM Structure and Return Adjacency Matrix
#'
#' @description
#' This function computes MDM scores and learns the structure of a dynamic Bayesian network,
#' returning its adjacency matrix using either `bnlearn::hc` or GOBNILP optimizer.
#'
#' @param data_input Data frame or array. Input multivariate time series (possibly with subjects).
#' @param method Character. "hc" to use hill-climbing from `bnlearn`, or "ipa" to use the GOBNILP optimizer.
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
  
  method <- match.arg(method, choices = c("hc", "ipa"))
  if(method == "hc"){
    if(verbose) message("Running hill-climbing with custom MDM score...")
    
    if(is.null(colnames(data_input))){
      colnames(data_input) <- paste0("V", seq_len(ncol(data_input)))
    }
    
    if(package_version(packageVersion("bnlearn"))>"4.9"){
      arg_score_name <- "custom-score"
    }else{
      arg_score_name <- "custom"
    }
    
    hill <- bnlearn::hc(
      x = data.frame(data_input),
      score = arg_score_name,
      fun = mdm_score_bn,
      args = list(nbf = nbf, method = "Brent", call = FALSE)
    )
    
    N <- ncol(data_input)
    var_names <- colnames(data_input)
    adj <- matrix(0, nrow = N, ncol = N, dimnames = list(var_names, var_names)
    )
    
    if(length(hill$arcs) > 0){
      for(z in seq_len(nrow(hill$arcs))){
        from <- hill$arcs[z, 1]
        to <- hill$arcs[z, 2]
        adj[from, to] <- 1
      }
    }
  } else if(method == "ipa"){
    if(is.null(gobnilp_path)) stop("You must provide gobnilp_path when using method = 'ipa'.")
    if(verbose) message("Computing MDM scores...")
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
