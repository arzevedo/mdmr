#' Delta selection
#'
#' @param dts The matrix with dataset; Number of timepoints X Number of nodes
#' @param m_ad Square Matrix Adjacent with dimension, Number of nodes  1 if edge exists; 0 otherwise
#' @param nbf the Log Predictive Likelihood will be calculate from this time point. It has to be a positive integer number; The default is 15
#' @param delta The vector with the sequence of all discount factors
#'
CDELT <- function(dts,m_ad,nbf=15,delta=seq(from=0.5, to=1.0, by=0.01)) {
      nd = length(delta)
      Nn = ncol(dts)
      Nt = nrow(dts)
      lpldet = array(0,dim=c(nd,Nn))

      for (k in 1:nd){
            for (i in 1:Nn){
                  # Initials:
                  p = sum(m_ad[,i])
                  if (m_ad[i,i] == 0) {p = p + 1}
                  Ft = array(1, dim=c(Nt,p))
                  aux = c(1:Nn)[m_ad[,i]>0]
                  aux2 = aux[aux!=i]
                  if (length(aux2)>0){ Ft[,2:(length(aux2)+1)] = dts[,aux2]}
                  Yt = dts[,i]
                  # DLM
                  a=dlm_filt(Yt, t(Ft), delta=delta[k])
                  lpldet[k,i]=sum(a$lpl[nbf:Nt])
            }
      }
      #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
      DF_hat = rep(0,Nn)
      for (i in 1:Nn){
            DF_hat[i] = na.omit(delta[lpldet[,i]==max(lpldet[,i],na.rm=TRUE)])[1]
      }
      result <- list(lpldet = lpldet, DF_hat = DF_hat)
      return(result)
}
