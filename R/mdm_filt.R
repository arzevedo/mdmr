#' Multiregression Dynamic Model (MDM) Filter
#'
#' @description
#' MDM with Filtering for unknown observational and state variances
#'
#' @param dts The matrix with dataset; Number of timepoints X Number of nodes
#' @param m_ad Square Matrix Adjacent with dimension = Number of nodes # 1 if edge exists; 0 otherwise
#' @param DF_hat Vector with delta that maximizes the LPL for each node with dimension = Number of nodes
mdm_filt <- function(dts, m_ad, DF_hat) {
      Nn = ncol(dts)
      Nt = nrow(dts)
      mt = vector(Nn, mode = "list")
      names(mt) <- colnames(m_ad)
      conections <- which(m_ad == 1, arr.ind = T)
      Ct = vector(Nn, mode = "list")
      Rt = vector(Nn, mode = "list")
      nt = vector(Nn, mode = "list")
      dt = vector(Nn, mode = "list")
      ft = vector(Nn, mode = "list")
      Qt = vector(Nn, mode = "list")
      ets = vector(Nn, mode = "list")
      lpl = vector(Nn, mode = "list")

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
            a=dlm_filt(Yt, t(Ft), delta=DF_hat[i])
            mt[[i]] = a$mt
            if(sum(m_ad[,i]) != 0){
                  rownames(mt[[i]]) <- c(
                        paste0("beta0_", colnames(m_ad)[i]),
                        paste0(
                              names(which(conections[,2] == i)),
                              "->", colnames(m_ad)[i]
                        )
                  )
            }

            Ct[[i]] = a$Ct
            Rt[[i]] = a$Rt
            nt[[i]] = a$nt
            dt[[i]] = a$dt
            ft[[i]] = a$ft
            Qt[[i]] = a$Qt
            ets[[i]] = a$ets
            lpl[[i]] = a$lpl
      }
      result <- list(mt=mt, Ct=Ct, Rt=Rt, nt=nt, dt=dt, ft=ft, Qt=Qt, ets=ets, lpl=lpl)
      return(result)
}
