#' Smooth dynamic linear models (DLMs)
#'
#'@description
#'DLM with Smoothing for unknown observational and state variances
#'
#' @param mt The matrix of posterior mean with dimension p X T
#' @param Ct The squared matrix of posterior variance with dimension p X p X T
#' @param Rt The squared matrix of prior variance with dimension p X p X T
#' @param nt The vector of prior hypermarameters of precision phi with length T
#' @param dt The vector of prior hypermarameters of precision phi with length T
#' @param Gt he matrix of state equation with dimension: p X p X T. The default is identity matrix block
#'
dlm_smoo <- function(mt, Ct, Rt, nt, dt, Gt = 0) {

      # defining objects
      if (is.vector(mt)){
            mt = array(mt, dim=c(1,length(mt)))
            Ct = array(Ct, dim=c(1,1,length(mt)))
            Rt = array(Rt, dim=c(1,1,length(Rt)))
      }
      if (Gt == 0){Gt = array(diag(nrow(mt)), dim=c(nrow(mt),nrow(mt),ncol(mt)))}
      p = nrow(mt) # the number of thetas
      Nt = ncol(mt) # the sample size
      smt = array(0, dim=c(p,Nt))
      sCt = array(0, dim=c(p,p,Nt))

      # in the last time point
      smt[,Nt] = mt[,Nt]
      sCt[,,Nt] = Ct[,,Nt]

      # for other time points
      for (i in (Nt-1):1){
            RSt = Rt[,,(i+1)]*nt[i]/dt[i]
            CSt = Ct[,,i]*nt[i]/dt[i]
            inv.sR = solvecov(RSt, cmax = 1e+10)$inv
            B = CSt %*% t(Gt[,,(i+1)]) %*% inv.sR
            smt[,i] = mt[, i] + B %*% (smt[,(i+1)] - Gt[,,(i+1)] %*% mt[,i])
            sCS = CSt + B %*% (sCt[,,(i+1)]*nt[Nt]/dt[Nt] - RSt) %*% t(B)
            sCt[,,i] = sCS * dt[Nt] / nt[Nt]
      }

      result <- list(smt=smt, sCt=sCt)
      return(result)
}
