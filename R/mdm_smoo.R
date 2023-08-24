#' Multiregression Dynamic Model (MDM) snooth
#'
#' @description
#' MDM with Smoothing for unknown observational and state variances
#'
#' @param mt List with the matrix of posterior mean with dimension p X T
#' @param Ct List with the squared matrix of posterior variance with dimension p X p X T
#' @param Rt List with the squared matrix of prior variance with dimension p X p X T
#' @param nt List with the vector of prior hypermarameters of precision phi with length T
#' @param dt List with the vector of prior hypermarameters of precision phi with length T
#'
mdm_smoo <- function(mt, Ct, Rt, nt, dt) {
      Nn = length(mt) # the number of nodes
      smt = vector(Nn, mode = "list")
      sCt = vector(Nn, mode = "list")
      SE = vector(Nn, mode = "list")

      for (i in 1:Nn){
            a=dlm_smoo(mt=mt[[i]], Ct=Ct[[i]], Rt=Rt[[i]], nt=nt[[i]], dt=dt[[i]])
            rownames(a$smt) <- rownames(mt[[i]])
            smt[[i]] = a$smt
            sCt[[i]] = a$sCt

            if(dim(sCt[[i]])[1] == 1){
                  SE[[i]] = qt(p = .975, df = nt[[i]])[length(qt(p = .975, df = nt[[i]]))] *
                        sqrt(a$sCt)
            } else {
                  aux = matrix(NA, ncol = dim(sCt[[i]])[1], nrow = dim(sCt[[i]])[3])
                  colnames(aux) <- paste0("SE_", rownames(mt[[i]]))
                  for(j in 1:dim(sCt[[i]])[1]) {

                        aux[, j] = qt(p = .975, df = nt[[i]])[length(qt(p = .975, df = nt[[i]]))] *
                              sqrt(sCt[[i]][j,j,])
                        SE[[i]] = aux
                  }
            }
      }
      result <- list(smt=smt, sCt=sCt, SE = SE)
      return(result)
}
