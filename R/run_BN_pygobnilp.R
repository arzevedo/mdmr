#' Network Structure Learning with pyGOBNILP
#'
#'@description
#'A function to learn structure of a Bayesian network (BN) from local scores
#'
#' @param score_file character File path to the local scores file
#' @param install_gob logical should the pyGOBNILP be install
#' @param palim_nodes numberical number of of nodes (variables).
#'
#' @import reticulate
run_BN_pygobnilp <- function(score_file = NULL, install_gob = FALSE,
                             palim_nodes = NULL){
      #library(reticulate)

      #if(require("reticulate")) stop("O pacote reticulate nao foi detectado.")

      # Criar um ambiente e instalar o pygobnilp
      #if(!("r-reticulate" %in% reticulate::virtualenv_list()))
      #      stop("Crie um ambiente virtual com a livraria pygobnilp")

      if(install_gob) reticulate::virtualenv_install("r-reticulate", "pygobnilp")

      # Definir o ambiente a ser utilizado
      reticulate::use_virtualenv("r-reticulate")

      # Importar o pacote e definir o modulo
      m <- reticulate::import("pygobnilp.gobnilp")$Gobnilp()

      # local_score run
      aux <- m$learn(local_scores_source = score_file, plot = FALSE,
                     palim = palim_nodes)


      dag_arrow <- unlist(reticulate::py_to_r(reticulate::tuple(m$arrow)))

      dag_arrow_char <- character()
      for(i in 1:length(dag_arrow)){
            dag_arrow_char[i] <- as.character(dag_arrow[[i]])
      }


      pattern <- "<gurobi\\.Var\\s(\\S+)->(\\S+)\\s\\(value\\s-?(\\d+\\.\\d+)\\)>"
      matches <- regmatches(dag_arrow_char, regexec(pattern, dag_arrow_char))

      # Create the data frame
      df <- data.frame(
            Father = sapply(matches, function(match) match[[2]]),
            Child = sapply(matches, function(match) match[[3]]),
            Value = as.numeric(sapply(matches, function(match) match[[4]]))
      )

      # adj_mat <- matrix(0, ncol = unique(df$Father), nrow = unique(df$Father))
      #
      # # Fill the matrix with values from the data frame
      # for (i in 1:nrow(adj_mat)) {
      #       adj_mat[df$Father[i], df$Child[i]] <- df$Value[i]
      # }

      # father <- unique(df$Father)
      # child <- unique(df$Child)
      #
      # adj_mat <- matrix(0, nrow = length(father), ncol = length(child))
      # for (i in 1:nrow(adj_mat)) {
      #       row_idx <- which(father == df$Father[i])
      #       col_idx <- which(child == df$Child[i])
      #       adj_mat[row_idx, col_idx] <- df$Value[i]
      # }
      #
      # rownames(adj_mat) <- father
      # colnames(adj_mat) <- child

      father <- unique(df$Father)
      child <- unique(df$Child)

      adj_mat <- matrix(NA, nrow = length(father), ncol = length(father),
                        dimnames = list(father, father))


      adj_mat[as.matrix(df[c(1,2)])] <- df$Value
      diag(adj_mat) <- 0

      result <- list(long = df, wide = adj_mat)

      return(result)

}
