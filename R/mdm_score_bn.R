#' Custom MDM Score Function for `bnlearn::hc`
#'
#' Calculates the Multiregression Dynamic Model (MDM) log predictive likelihood score 
#' to be used as a custom scoring function in `bnlearn::hc()`.
#'
#' @param node A character string representing the target node.
#' @param parents A character vector of parent node names.
#' @param data A data frame of the observed multivariate time series.
#' @param args A list of arguments including:
#'   \itemize{
#'     \item \code{nbf}: Numeric, the burn-in time point.
#'     \item \code{method}: Optimization method for delta (e.g., "Brent").
#'     \item \code{call}: Logical, whether to transform the data internally.
#'   }
#'
#' @return The log predictive likelihood score (numeric).
mdm_score_bn <- function(node, parents, data, args) {
  n_n = dim(data)[2]
  m_ad = array(0, dim=c(n_n,n_n))
  dimnames(m_ad) = list(colnames(data),colnames(data))
  #m_ad = m_ad[dimnames(m_ad)[[1]][order(nchar(dimnames(m_ad)[[1]]), dimnames(m_ad)[[1]])],dimnames(m_ad)[[2]][order(nchar(dimnames(m_ad)[[2]]), dimnames(m_ad)[[2]])]]
  m_ad[parents,node] = 1
  nbf= args$nbf
  method = args$method
  call = args$call
  #nd = length(delta)
  Nn = ncol(data)
  Nt = nrow(data)
  i = grep(node, colnames(m_ad))[1]
  #lpldet = rep(0,Nn)
  if (call){
    data_input = data
    data = array(0, dim=c(Nt, Nn))
    for(i in 1:Nn){
      data[,i] <- data_input[,i]
    }
  }
  
  res = list()
  k<-optim(logpl,par=0.5,dts=data,m_ad=m_ad,i=i,method=method,lower=0,upper=1)
  res<-k
  lpldet=-res$value
  
  
  #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
  DF_hat = rep(0,Nn)
  for (i in 1:Nn){
    DF_hat[i] = res$par
  }
  #result <- list(lpldet = lpldet, DF_hat = DF_hat)
  #return(result)
  return(lpldet)
}