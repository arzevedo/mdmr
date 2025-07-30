#' Local Predictive Log-Likelihood Function
#'
#' Computes the negative log predictive likelihood used in model scoring 
#' for a given node in a dynamic linear model (DLM), conditioned on its parents.
#'
#' @param dts A matrix of time series data (observations × nodes).
#' @param m_ad A binary adjacency matrix indicating parent relationships.
#' @param delta Discount factor used in the DLM.
#' @param i Index of the target node.
#'
#' @return A single numeric value representing the negative log predictive likelihood.
#' @keywords internal
logpl<-function(dts,m_ad,delta,i){
  Nt = nrow(dts)
  Nn = ncol(dts)
  nbf=15
  p = sum(m_ad[,i])
  if (delta>1|delta<0){
    return(Inf)
  }
  if (m_ad[i,i] == 0) {p = p + 1}
  Ft = array(1, dim=c(Nt,p))
  aux = c(1:Nn)[m_ad[,i]>0]
  aux2 = aux[aux!=i]
  if (length(aux2)>0){
    if (typeof(dts[,aux2])=='list'){
      Ft_aux = matrix(unlist(dts[,aux2]),nrow=Nt,ncol=length(aux2))
      Ft[,2:(length(aux2)+1)] = Ft_aux
    }else{
      Ft[,2:(length(aux2)+1)] = dts[,aux2]
    }
  }
  Yt = dts[,i]
  # DLM
  a=dlm_filt(Yt, t(Ft), delta=delta)
  lpldet=sum(a$lpl[nbf:Nt])
  return(-lpldet)
}