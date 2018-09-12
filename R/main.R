#' Generalized fiducial inference
#'
#' @param dat data
#' @param FE design matrix of fixed effects
#' @param RE random effects
#' @param N desired number of simulations
#' @param thresh threshold
#'
#' @return A list with two elements:
#' 
#' @export
#' 
#' @importFrom stats rnorm rchisq
#' @importFrom Matrix rankMatrix
#'
#' @examples
#' h <- 0.01
#' dat <- data.frame(npk$yield-h, npk$yield+h)
#' FE <- as.matrix(rep(1,nrow(dat))) # intercept
#' RE <- data.frame(block = npk$block)
#' gfi <- gfimm(dat, FE, RE, N=500)
gfimm <- function(dat, FE, RE, N, thresh=N/2){
  n <- nrow(dat)
  #random_design  <-  1 
  L <- dat[,1L] 
  U <- dat[,2L] 
  fe <- ncol(FE) 
  break_point <- 10L
  ########--------SET-UP RANDOM EFFECTS
  # RE is declared.
  re <- ncol(RE)+1L 
  E <- integer(re) # E[i] = number of levels of i-th random effect
  E[re] <- n 
  for(i in 1L:(re-1L)){
    E[i] <- nlevels(RE[,i]) 
  }
  ESUM <- cumsum(E)
  RE2 <- cbind(RE, factor(1L:n)) #Adds the error effect 
  RE <- NULL 
  for(i in 1L:re){ # Builds an indicator RE matrix for the effects
    re_levels <- levels(RE2[,i])
    for(j in 1L:E[i]){
      temp1 <- which(RE2[,i]==re_levels[j]) 
      temp2 <- integer(n) 
      temp2[temp1] <- 1L 
      RE <- cbind(RE,temp2)
    }
  } 
  
  RE2 <- sapply(RE2, as.integer)
  
  Dim <- fe+re  # Dimension of the space
  
  ########--------SET-UP ALGORITHM OBJECTS
  Z <- weight <- vector(mode = "list", length = re) # Particles and Weights
  ESS <- rep(N,n) # Effective sample size
  VC <- integer(N) # Number of vertices
  CC <- vector(mode = "list", length = N) # constraints 
  C <- NULL # initial constraints 
  VT <- vector("list", N) # vertices
  
  ########--------SAMPLE ALL Z's/SET-UP WEIGHTS
  A <- array(NA_real_, dim = c(n, re, N)) 
  for(i in 1L:re){  #Sample all the Z's and set-up the weights
    Z[[i]] <- matrix(rnorm(E[i]*N), ncol=N) 
    weight[[i]] <- matrix(1, nrow=E[i], ncol=N)
    M <- RE[, (ESUM[i]-E[i]+1L):ESUM[i]] %*% Z[[i]]
    for(j in 1L:N){
      A[,i,j] <- M[,j]
    } 
  }      
  
  AA <- cbind(FE, A[,,1L])	
  AT <- NULL 
  r <- 0L   
  if(rankMatrix(AA) != Dim){
    stop("Design is not of full rank")
  }else{
    for(k in 1L:n){
      A_temp <- rbind(AT, AA[k,])  
      if(rankMatrix(A_temp) > r){
        AT <- A_temp
        r <- rankMatrix(AT)
        C <- c(C, k) 
      }
    }
  }
  K <- c(1L:n)[!is.element(1L:n, C)] 
  K_start <- C
  Z[[re]][K,] <- 0  #remove error particles not selected in initialization
  # K will be the values left to be sampled
  #} # C=K_start are the indexes of Z[[residual]] we keep 
  
  ########--------FIND INITIAL VERTICES
  # input : C=K_start
  # one line of USE = one combination of Dim indices such that the corresponding rows of [[FE A] ; [-FE -A]] are free
  USE <- as.matrix(expand.grid(lapply(1L:Dim, function(i) c(0L,n) + C[i]), 
                               KEEP.OUT.ATTRS=FALSE))
  # now calculate vertices of Q1 in V 
  b <- c(U,-L) 
  twoPowerDim <- as.integer(2^Dim) # number of vertices
  FEFE <- rbind(FE,-FE)
  for(j in 1L:N){
    V <- matrix(NA_real_, nrow=Dim, ncol=twoPowerDim) 
    Aj <- A[,,j]
    AA <- cbind(FEFE, rbind(Aj, -Aj))
    for(ii in 1L:twoPowerDim){ 
      II <- USE[ii,]  # constraints to use are by row
      V[,ii] <- solve(AA[II,],b[II])  # selection in Q ; the V function
    }
    VT[[j]] <- V 
  }
  CC <- rep(list(t(USE)), N) #constraints are the same for all N particles
  VC <- rep(twoPowerDim, N) # this is (ncol(CC[[1]]), ..., ncol(CC[[N]])) 
  
  ########--------MAIN ALGORITHM
  K_n <- as.integer(ceiling(length(K)/break_point)) # split K for sampling 
  K_temp <- vector("list", K_n)
  for(i in 1L:(K_n-1L)){
    K_temp[[i]] <- K[((i-1L)*break_point+1L):(i*break_point)]
  }
  K_temp[[K_n]] <- K[((K_n-1L)*break_point+1L):length(K)] 
  K1 <- NULL
  for(k_n in 1L:K_n){
    K1 <- c(K1, K_temp[[k_n]])
    for(k in K_temp[[k_n]]){
      if(k_n > 1L){
        for(i in 1L:re){
          Zi <- Z[[i]]
          Z[[i]] <- rbind(Zi, matrix(rnorm(N*(E[i]-nrow(Zi))),ncol=N))  
        }
      }
      for(i in 1L:N){
        VT1 <- VTi <- VT[[i]] # vertex values
        VT2 <- VT1[Dim,]
        VT1 <- VT1[-Dim,] 
        if(fe>0){
          Z1 <- FE[k,]
        }else{
          Z1 <- NULL
        }
        for(j in 1L:re){
          # if(RE2[k,j]==0L){ # RE2 never zero!
          #   Z1 <- c(Z1,0)
          # }else{ 
          Z1 <- c(Z1, Z[[j]][RE2[k,j],i])
          # }
        }
        Z1 <- Z1[-Dim] # remove column to be sampled
        VTsum <- c(Z1 %*% VT1)
        fidsample <- fid_sample(VT2, VTsum, U[k], L[k]) ###Sample	
        ZZ <- fidsample["ZZ"]
        wt <- fidsample["wt"]
        Z[[re]][k,i] <- ZZ
        weight[[re]][k,i] <- wt
        VTsum <- VTsum + ZZ*VT2
        VT1 <- VTi 
        CC1 <- CC[[i]]
        fidvertex <- fid_vertex(VT1, CC1, VTsum, U, L, Dim, k, n)
        VC[i] <- fidvertex$vert # this is ncol(CCtemp)  
        CC[[i]] <- fidvertex$CCtemp
        # if(ncol(CC[[i]]) != VC[i]){ #check
        #   cat("ncol(fidvertex$CCtemp) != fidvertex$vert", "\n")
        # }
        VT[[i]] <- fidvertex$VTtemp
        if(fidvertex$vert == 0L){ #check
          cat("fidvertex$vert==0L", "\n") # never occurs in the examples I've tried
          weight[[re]][k,i] <- 0
        }
      }
      WT <- apply(weight[[re]], 2L, cumprod)  #only last re is restricted
      WT <- WT[nrow(WT),]
      if(sum(WT)==0){
        stop("Error: possible underflow")
      }
      WT <- WT/sum(WT)
      ESS[k] <- 1/c(crossprod(WT))
      #---------------------------------------------------------------Resample
      if(ESS[k]<thresh && k<max(K)){
        N_sons <- integer(N)  
        # generate the cumulative distribution
        dist <- cumsum(WT) 
        aux <- runif(1L) # sample uniform rv in [0 1]
        u <- (aux + 0L:(N-1L))/N # or is it  aux + (0L:(N-1L))/N ? no
        j <- 1L 
        for(i in 1L:N){
          while(u[i]>dist[j]){ 
            j <- j+1L 
          }
          N_sons[j] <- N_sons[j]+1L 
        }
        JJ <- unique(c(1L:k, C)) # K_start and integers <=k
        II <- c(1L:n)[-JJ]
        ZZ <- vector("list", re)			
        VCVC <- integer(N) 
        CCCC <- VTVT <- vector("list", N)
        for(i in 1L:N){  
          if(N_sons[i]){
            VCtemp <- rep(VC[i], N_sons[i]) 
            Ztemp <- vector("list", re)
            VTtemp <- vector("list", N_sons[i])
            copy <- N_sons[i]-1L  # to be resampled
            VTi <- VT[[i]]
            for(ii in 1L:N_sons[i]){
              VTtemp[[ii]] <- VTi # copy original vertices
              # could do VTtemp[1L:N_sons[i]] <- VTi ? rather <- rep(list(VTi), N_sons[i])
            }
            for(ii in 1L:re){
              Ztemp[[ii]] <- Z[[ii]][,rep(i,N_sons[i])] #copy Z
            }			
            if(copy){
              ord <- sample.int(re)  #Order to resample. Each re will be resampled.
              for(kk in ord){
                for(ii in 1L:copy){
                  XX <- NULL
                  for(jj in 1L:re){
                    XX <- cbind(XX, 
                                RE[,(ESUM[jj]-E[jj]+1L):ESUM[jj]] %*% 
                                  Ztemp[[jj]][,ii]) 
                  }
                  XX <- XX[JJ,-kk] #remove column of effect to be resampled
                  ### temp <- which(RE2[JJ,kk]!=0L)  #find which levels of kk have been sampled          
                  ### RE2 is never 0 !
                  temp <- seq_along(JJ)
                  #
                  Z1 <- Ztemp[[kk]][unique(RE2[JJ[temp],kk]),ii] #Z being resampled
                  CO2 <- RE[JJ,(ESUM[kk]-E[kk]+1L):ESUM[kk]] 
                  CO2 <- CO2[,which(colSums(CO2)!=0L)] #levels not sampled yet
                  # Z1==0 <=> Z's not sampled
                  Z00 <- which(Z1 != 0) #These are the levels with Z for effect kk
                  Z1 <- Z1[Z00]
                  if(fe>0){
                    XX <- cbind(FE[JJ,],XX) 
                  }
                  MAT <- cbind(-XX,CO2)
                  kernel <- nullSpace(MAT)
                  # qrMAT <- qr(t(MAT), tol=1e-10)
                  # if(qrMAT$rank < ncol(MAT)){ # faire qr pour avoir le rank, et utilise le r?sultat pour le null space (voir le code de Null)
                  # NUL <- qr.Q(qrMAT, complete = TRUE)[, -seq_len(qrMAT$rank), drop=FALSE]
                  if(kernel$rank < ncol(MAT)){
                    NUL <- kernel$kernel
                    n1 <- NUL[1L:(nrow(NUL)-ncol(CO2)), , drop=FALSE] # eta1 in [HC]  
                    n2 <- NUL[(nrow(NUL)-ncol(CO2)+1L):nrow(NUL), , drop=FALSE] # eta2 in [HC], up to normalization 
                    # GS <- gramSchmidt(n2)
                    # O2 <- GS$Q
                    QR <- qr(n2)
                    O2 <- qr.Q(QR)
                    R <- qr.R(QR)
                    O1 <- tsolveAndMultiply(R, n1) # n1 %*% backsolve(R, diag(nrow(R)))
                    a <- t(O2) %*% Z1 # C in [HC] but with n2 and not O2
                    tau <- Z1 - c(O2 %*% a)
                    b <- sqrt(crossprod(tau)[1L])
                    tau <- tau/b
                    rank.O2 <- ncol(O2) 
                    bb <- sqrt(rchisq(1L,length(Z1)-rank.O2)) 
                    bbb <- b/bb	  					
                    aa <- rnorm(rank.O2) 
                    MM3 <- O2 %*% aa + bb*tau
                    Ztemp[[kk]][Z00,ii] <- c(MM3) 
                    vert <- c(1L:Dim)[-(fe+kk)]
                    # M1 <- CO2 %*% Z1
                    M2 <- c(O1 %*% (bbb*aa-c(a)))
                    # M3 <- CO2 %*% MM3
                    VTtemp_ii <- VTtemp[[ii]]
                    for(jj in 1L:VC[i]){ 
                      # check1 <- XX %*% VTtemp[[ii]][vert,jj] + 
                      #   VTtemp[[ii]][fe+kk,jj] * M1 
                      VTtemp[[ii]][vert,jj] <- VTtemp_ii[vert,jj] - 
                        VTtemp_ii[fe+kk,jj] * M2 
                      VTtemp[[ii]][fe+kk,jj] <- VTtemp_ii[fe+kk,jj]*bbb
                      # check2 <- XX %*% VTtemp[[ii]][vert,jj] + 
                      #   VTtemp[[ii]][fe+kk,jj] * M3
                      # if(!all.equal(check1, check2)){
                      #   cat("check1!=check2", "\n")
                      # }
                    }
                  }else{
                    b <- sqrt(crossprod(Z1)[1L]) 
                    tau <- Z1/b 
                    bb <- sqrt(rchisq(1L,length(Z1))) 
                    Ztemp[[kk]][Z00,ii] <- bb*tau 
                    vert <- c(1L:Dim)[-(fe+kk)] 
                    # VTtemp_ii <- VTtemp[[ii]]
                    # for(jj in 1L:VC[i]){
                    #   VTtemp[[ii]][fe+kk,jj] <- VTtemp_ii[fe+kk,jj]*b/bb 
                    # }
                    VTtemp[[ii]][fe+kk,1L:VC[i]] <- 
                      VTtemp[[ii]][fe+kk,1L:VC[i]]*b/bb
                  }
                }
              }
              #}
            } 
            for(ii in 1L:re){
              ZZ[[ii]] <- cbind(ZZ[[ii]], Ztemp[[ii]])
            }
            d <- sum(N_sons[seq_len(i-1L)])
            VCVC[(d+1L):sum(N_sons[1L:i])] <- VCtemp 
            for(kk in 1L:N_sons[i]){ 
              VTVT[[kk+d]] <- VTtemp[[kk]]
              CCCC[[kk+d]] <- CC[[i]]
            }
          }
        }
        Z <- ZZ
        VT <- VTVT
        VC <- VCVC
        CC <- CCCC
        weight[[re]] <- matrix(1, nrow=E[re], ncol=N) #assign weights of error matrix to 1
      }
    } # ends resampling for k=K1
    #----------------------------------------------------determine signs
    signs <- matrix(0L, nrow=re, ncol=N)    
    for(i in 1L:N){
      VTi <- VT[[i]]
      for(j in 1L:re){
        if(all(VTi[fe+j,]>0)){ #i.e. all are positive
          signs[j,i] <- 1L
        }
        else if(all(VTi[fe+j,]<0)){ #i.e. all are negative
          signs[j,i] <- -1L 
        }
      }
    }
    #----------------------------------------------------FINAL RESAMPLE			
    ZZ <- nn <- vector("list", re)
    VTVT <- vector("list", N)
    n1 <- sort(c(K1, K_start)) # = 1:n for last k_n 
    for(ii in 1L:re){
      nn[[ii]] <- unique(RE2[n1,ii])  
    }
    lengths.nn <- lengths(nn, use.names = FALSE)
    for(i in 1L:N){
      Ztemp <- vector("list", re)
      VTtemp <- VT[[i]] 
      for(ii in 1:re){
        Ztemp[[ii]] <- Z[[ii]][nn[[ii]],i]  #copy Z
      }
      ord <- sample.int(re)
      for(kk in ord){
        XX <- NULL
        eff <- c(1L:re)[-kk] 
        for(jj in eff){
          XX <- cbind(
            XX, 
            RE[,(ESUM[jj]-E[jj]+1L):(ESUM[jj]-E[jj]+lengths.nn[jj])] %*% 
              Ztemp[[jj]]
          ) 
        }
        #temp <- which(RE2[,kk]!=0)  #find which levels of kk have been sampled
        # RE2 never 0 ! 
        Z1 <- Ztemp[[kk]]  #Z being resampled
        CO2 <- RE[,(ESUM[kk]-E[kk]+1L):(ESUM[kk]-E[kk]+lengths.nn[kk])] 
        ## Z1==0 <=> Z's not sampled
        # Z00 <- which(Z1 != 0) #These are the levels with Z for effect kk
        # if(length(Z00) < length(Z1)){
        #   cat("resampling - some Z1's not 0") # this never occurs
        # }
        # Z1 <- Z1[Z00]
        Z00 <- seq_along(Z1)
        if(fe>0){
          XX <- cbind(FE, XX)
        }
        MAT <- cbind(-XX, CO2)
        kernel <- nullSpace(MAT)
        # qrMAT <- qr(t(MAT), tol=1e-10)
        # if(qrMAT$rank < ncol(MAT)){ # faire qr pour avoir le rank, et utilise le r?sultat pour le null space (voir le code de Null)
        # NUL <- qr.Q(qrMAT, complete = TRUE)[, -seq_len(qrMAT$rank), drop=FALSE]
        if(kernel$rank < ncol(MAT)){
          NUL <- kernel$kernel
          n1 <- NUL[1L:(nrow(NUL)-ncol(CO2)), , drop=FALSE]
          n2 <- NUL[(nrow(NUL)-ncol(CO2)+1L):nrow(NUL), , drop=FALSE]
          # GS <- gramSchmidt(n2)
          # O2 <- GS$Q
          # O1 <- n1 %*% backsolve(GS$R, diag(nrow(GS$R)))
          QR <- qr(n2)
          O2 <- qr.Q(QR)
          R <- qr.R(QR)
          O1 <- tsolveAndMultiply(R, n1) # n1 %*% backsolve(R, diag(nrow(R)))
          a <- t(O2) %*% Z1
          tau <- Z1 - c(O2%*%a)
          b <- sqrt(crossprod(tau)[1L])
          tau <- tau/b
          rank.O2 <- ncol(O2) 
          bb <- sqrt(rchisq(1L,length(Z1)-rank.O2))
          bbb <- b/bb	 
          aa <- rnorm(rank.O2)
          Ztemp[[kk]][Z00] <- O2%*%aa + bb*tau # Z00 useless
          vert <- c(1L:Dim)[-(fe+kk)]
          M2 <- O1 %*% (bbb*aa-c(a))
          for(jj in 1L:VC[i]){
            VTtemp[vert,jj] <- VTtemp[vert,jj]- VTtemp[fe+kk,jj] * M2 
            VTtemp[fe+kk,jj] <- VTtemp[fe+kk,jj] * bbb 
          }
        }else{
          b <- sqrt(crossprod(Z1)[1L])
          tau <- Z1/b
          bb <- sqrt(rchisq(1L,length(Z1)))
          Ztemp[[kk]][Z00] <- bb*tau # Z00 useless
          vert <- c(1L:Dim)[-(fe+kk)]
          bbb <- b/bb
          # for(jj in 1L:VC[i]){
          #   VTtemp[fe+kk,jj] <- VTtemp[fe+kk,jj]*bbb
          # }
          VTtemp[fe+kk,1L:VC[i]] <- VTtemp[fe+kk,1L:VC[i]]*bbb
        }
      }    
      for(ii in 1L:re){
        ZZ[[ii]] <- cbind(ZZ[[ii]], Ztemp[[ii]]) 
      }
      VTVT[[i]] <- VTtemp  
    }
    Z <- ZZ 
    VT <- VTVT
    #----------------------------------------------------flip negatives			
    for(i in 1L:N){
      VTtemp.i <- VT[[i]]
      for(j in 1L:re){ #only look at random effects
        if(signs[j,i] == -1L){
          VTtemp.i[fe+j,] <- -VTtemp.i[fe+j,] #only need to flip the negatives
          Z[[j]][,i] <- -Z[[j]][,i]             
        }
      }
      VT[[i]] <- VTtemp.i
    }    
    if(k_n == K_n){ #if finished pick coordinates			
      #pick the coordinates
      unif <- matrix(runif(Dim*N), nrow=Dim, ncol=N) 
      # VT_end <- matrix(0, nrow=Dim, ncol=N) 
      # for(i in 1L:N){
      #   VTi <- VT[[i]]
      #   for(j in 1L:Dim){
      #     if(unif[j,i] <= 0.5){
      #       if(j <= fe){
      #         VT_end[j,i] <- min(VTi[j,])
      #       }else{
      #         VT_end[j,i] <- max(min(VTi[j,]),0) 
      #       }
      #     }else{
      #       if(j <= fe){
      #         VT_end[j,i] <- max(VTi[j,])
      #       }else{
      #         VT_end[j,i] <- max(max(VTi[j,]),0) 
      #       }
      #     }
      #   }
      # }
      VT_end <- pickCoordinates(Dim, N, fe, VT, unif)
      VERTEX <- VT_end 
      WEIGHT <- WT 
    }
  }
  return(list(VERTEX=VERTEX,WEIGHT=WEIGHT))
}
