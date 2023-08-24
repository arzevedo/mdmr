#' Multiregression Dynamic Model (MDM) Score
#'
#' @description
#' MDM scores
#'
#' @param data_input The DataFrame that holds the multivariate time series data
#' @param nbf Numeric, The Log Predictive Likelihood will be calculate from this time point. It has to be a positive integer number; The default is 15
#' @param delta Vector, Sequence of all discount factors. The default is seq(from=0.5, to=1.0, by=0.01)
#' @param GOLB_print Logical, if TRUE, save a file that can be used in the GOBNILP optimizer
#' @param subjects_length Numeric, how many subjects are contained in dataset
#'
mdm_score <- function(data_input, nbf=15, delta=seq(from=0.5, to=1.0, by=0.01),
                      GOLB_print = FALSE, subjects_length = 1) {

      Nt = dim(data_input)[1] #number of time points
      Nn = dim(data_input)[2] #number of nodes
      Ns = subjects_length #number of subjects

      dts = array(0, dim=c(Ns,Nt,Nn))

      if(subjects_length == 1){
            for(i in 1:Nn){
                  dts[,,i] <- data_input[,i]
            }
      } else {
            for(j in 1:Ns){
                  for(k in 1:Nn){
                        dts[j,,k] <- data_input[,k]
                  }
            }
      }
      # Add names to nodes
      dimnames(dts) <- list(1:Ns, 1:Nt, colnames(data_input))

      Ns = dim(dts)[1] # the number of subjects
      Nt = dim(dts)[2] # the number of timepoints
      Nn = dim(dts)[3] # the number of nodes
      Nd =  Nn*2^(Nn-1)# the number of possible parent-child

      delt_hat = array(0,dim=c(Ns,Nd)) # the discount factor chosen for each possibility
      lpl = array(0,dim=c(Ns,Nd)) # scores
      par_chil = array(-9,dim=c(Ns,Nd,(Nn+2))) #col1=subject; col2=model; col3={child,number of parents,parents}

      # generating all possible combinations
      cc = vector(Nn, mode = "list")
      for (i in 1:Nn){
            cc[[i]] = combn(c(1:Nn),i)
      }

      # finding the scores ----

      for (s in 1:Ns){
            #cat("loop subject ",s,"/",Ns, "\n")
            Np=0
            for (i in seq(from=1,to=Nn,by=2)){
                  if(subjects_length > 1){
                  cat("looping subject ",s,"/",Ns," and ",i,"/",Nn,"\n")
                  } else {
                        cat("looping node ",i,"/",Nn,"\n")
                  }
                  for (j in 1:ncol(cc[[i]])){
                        # Matrix Adjacent: 1 for edge; 0 otherwise
                        m_ad = diag(Nn)
                        p=cc[[i]][,j]
                        m_ad[p,]=1
                        # choosing the Discount Factor and finding the scores
                        delt_dag = CDELT(dts[s,,],m_ad, nbf=nbf, delta=delta)
                        delt_hat[s,(Np+1):(Np+Nn)]=delt_dag$DF_hat
                        #aux=t(delt_dag$lpldet)
                        #lpl[s,(Np+1):(Np+Nn)]=diag(aux[,max.col(aux)])
                        lpl[s,(Np+1):(Np+Nn)]=apply(delt_dag$lpldet,2,max,na.rm=TRUE) #deleting the NA from lpl
                        # child
                        par_chil[s,(Np+1):(Np+Nn),1] = seq(1:Nn)-1
                        # number of parents
                        par_chil[s,(Np+1):(Np+Nn),2] = i
                        par_chil[s,c((Np+1):(Np+Nn))[p],2] = i-1
                        # parents
                        par_chil[s,(Np+1):(Np+Nn),3:(2+i)] = t(array((p-1), dim=c(i,Nn)))
                        if (i==1){par_chil[s,c((Np+1):(Np+Nn))[p],3:(2+i)] = -9}
                        else {diag(par_chil[s,c((Np+1):(Np+Nn))[p],3:(2+i)]) = -9}

                        Np=Np+Nn
                  }
            }
      }

      DF_hat= array(-9,dim=c(Ns,Nd,(Nn+2)))
      DF_hat[,,1]=delt_hat
      DF_hat[,,2]=par_chil[,,1]
      DF_hat[,,3:(Nn+2)]=par_chil[,,3:(Nn+2)]
      dimnames(DF_hat)=list(c(1:Ns),c(1:Nd),c("DF","node",1:Nn))

      # creating the file with structure of James'programm
      all.score = list()
      for (s in 1:Ns){
            all.score[[s]] = array(" ",dim=c((Nd+Nn+1),(Nn+1)))
            all.score[[s]][1,1] = Nn
            a = cbind(lpl[s,],par_chil[s,,])
            b=a[order(a[,2]),]
            more=2
            j=0
            for (i in 1:Nd){
                  if (j==b[i,2]){
                        all.score[[s]][more,1:2] = c(j,Nd/Nn)
                        j=j+1
                        more=more+1
                  }
                  d = b[i,c(1,3:(Nn+3))]
                  all.score[[s]][more,1:length(d[d!=-9])] = d[d!=-9]
                  more = more + 1
            }
      }

      all.score_dic <- all.score[[1]]
      for(i in 1:Nn){
            all.score_dic[all.score_dic == i-1] <- gsub(x = dimnames(dts)[[3]][i],
                                                        " ", replacement =  "")
      }

      if(GOLB_print){
            write.table(x = all.score_dic,
                        file = paste0("mdm_score_", format(Sys.time(), "%d_%b_%Y")),
                        quote = FALSE, row.names = FALSE,
                        col.names = FALSE)

      }

      result <- list(all.score=all.score_dic, DF_hat=DF_hat)
      return(result)
}
