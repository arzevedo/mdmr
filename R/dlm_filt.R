#' Dynamic linear models (DLMs) filter
#'
#'@description
#'DLM with Filtering for unknown observational and state variances
#'
#' @param Yt The vector of observed time series with length T
#' @param Ft The matrix of covariates with dimension: number of thetas (p) X sample size (T)
#' @param delta Discount factor | Wt=Ctx(1-delta)/delta
#' @param Gt The matrix of state equation with dimension: p X p X T. The default is identity matrix block
#' @param m0 The vector of prior mean at time t=0 with length p. The default is non-informative prior, with zero mean.
#' @param CS0 The squared matrix of prior variance - C*0 | C*0Vt = C0, with length p. The default is non-informative prior, with prior variance equal to 3 times the observed variance.
#' @param n0 The prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative prior, with value of 0.001. n0 has to be higher than 0.
#' @param d0 The prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative prior, with value of 0.001. n0 has to be higher than 0.
#' @export
dlm_filt <- function(Yt, Ft, delta, Gt = array(diag(nrow(Ft)), dim=c(nrow(Ft),nrow(Ft),length(Yt))), m0 = rep(0,nrow(Ft)), CS0 = 3*diag(nrow(Ft)), n0 = 0.001, d0 = 0.001) {

      # defining objects
      p = nrow(Ft) # the number of thetas
      Nt = length(Yt)+1 # the sample size + t=0
      if (n0 == 0){
            n0 = 0.001
            warning("n0 is 0.001")
      }
      Y = rep(0, Nt)
      Y[2:Nt] = Yt
      F = array(0, dim=c(p,Nt))
      F[,2:Nt] = Ft
      G = array(0, dim=c(p,p,Nt))
      G[,,2:Nt] = Gt
      mt = array(0, dim=c(p,Nt))
      mt[,1] = m0
      Ct = array(0, dim=c(p,p,Nt))
      Ct[,,1] = CS0*d0/n0
      Rt = array(0, dim=c(p,p,Nt))
      nt = rep(0, Nt)
      nt[1] = n0
      dt = rep(0, Nt)
      dt[1] = d0
      ft = rep(0, Nt)
      Qt = rep(0, Nt)
      ets = rep(0, Nt)
      lpl = rep(0, Nt)

      for (i in 2:Nt){

            # Posterior at {t-1}: (theta_{t-1}|y_{t-1}) ~ t_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1}xd_{t-1}/n_{t-1}]
            # Prior at {t}: (theta_{t}|y_{t-1}) ~ t_{n_{t-1}}[a_{t}, R_{t}]
            at = G[,,i] %*% mt[,(i-1)]
            RSt = G[,,i] %*% (Ct[,,(i-1)]*nt[(i-1)]/dt[(i-1)]) %*% t(G[,,i]) / delta
            Rt[,,i] = RSt * dt[(i-1)] / nt[(i-1)]

            # One-step forecast: (Y_{t}|y_{t-1}) ~ t_{n_{t-1}}[f_{t}, Q_{t}]
            ft[i] = t(F[,i]) %*% at
            QSt = t(F[,i]) %*% RSt %*% F[,i] + 1
            Qt[i] = QSt * dt[(i-1)] / nt[(i-1)]
            et = Y[i] - ft[i]
            ets[i] = et / sqrt(Qt[i])

            # Posterior at t: (theta_{t}|y_{t}) ~ t_{n_{t}}[m_{t}, C_{t}]
            At = Rt[,,i] %*% F[,i] / Qt[i]
            mt[,i] = at + At * et
            nt[i] = nt[(i-1)] + 1
            dt[i] = dt[(i-1)] + (et^2) / QSt
            CSt = RSt - (At %*% t(At)) * QSt[1]
            Ct[,,i] = CSt * dt[i] / nt[i]

            # Log Predictive Likelihood
            lpl[i] <- lgamma((nt[i-1]+1)/2)-lgamma(nt[i-1]/2)-0.5*log(pi*nt[i-1]*Qt[i])-((nt[i-1]+1)/2)*log(1+(1/nt[i-1])*et^2/Qt[i])

      }

      result <- list(mt=mt[,2:Nt], Ct=Ct[,,2:Nt], Rt=Rt[,,2:Nt], nt=nt[2:Nt], dt=dt[2:Nt], ft=ft[2:Nt], Qt=Qt[2:Nt], ets=ets[2:Nt], lpl=lpl[2:Nt])

      return(result)
}
