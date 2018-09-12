fid_vertex <- function(VT1, CC1, VTsum, U, L, Dim, k, n){
  p <- length(VTsum) 
  whichl <- VTsum >= L[k] #vertices that satisfy lower constraint
  whichu <- VTsum <= U[k] #vertices that satisfy upper constraint
  both <- whichl & whichu
  checkl <- !whichl
  checku <- !whichu
  CCtemp <- matrix(integer(0L), nrow=Dim, ncol=0L) 
  VTtemp <- matrix(numeric(0L), nrow=Dim, ncol=0L)
  vert <- 0L 
  CA <- CC1[, checkl, drop=FALSE]
  CB <- CC1[, whichl, drop=FALSE] 
  if(any(checkl)){ #i.e. they do not all satisfy the lower constraint
    l <- sum(checkl) 
    if(l < p){
      INT <- matrix(0L, nrow=2L*n, ncol=l)
      for(ll in 1L:l){ #check lower constraints first
        INT[CA[,ll],ll] <- 1L 
      }
      VTsum_l <- VTsum[whichl]
      VTsum_ll <- VTsum[checkl]
      VT1_l <- VT1[, checkl, drop=FALSE]
      VT1_ll <- VT1[, whichl, drop=FALSE]
      for(ii in seq_len(p-l)){
        INT2 <- INT[CB[,ii], , drop=FALSE] 
        use <- which(colSums(INT2) == Dim-1L)
        for(dd in use){
          inter <- CB[INT2[,dd]==1L,ii] # this will be intersection
          vert <- vert + 1L
          CCtemp <- cbind(CCtemp, c(inter,k+n)) # need to add n indicating lower constraint
          lambda <- (L[k]-VTsum_l[ii])/(VTsum_ll[dd]-VTsum_l[ii])
          VTtemp <- cbind(VTtemp,
                          lambda*VT1_l[,dd] + (1-lambda)*VT1_ll[,ii])
        }
      }
    }
  }
  CA <- CC1[, checku, drop=FALSE] 
  CB <- CC1[, whichu, drop=FALSE]
  if(any(checku)){ # i.e. they do not all satisfy the lower constraint
    l <- sum(checku)
    if(l < p){
      INT <- matrix(0L, nrow=2L*n, ncol=l) 
      for(ll in 1L:l){ # check lower constraints first
        INT[CA[,ll],ll] <- 1L 
      }
      VTsum_u <- VTsum[whichu]
      VTsum_uu <- VTsum[checku]
      VT1_u <- VT1[,checku,drop=FALSE]
      VT1_uu <- VT1[,whichu,drop=FALSE]
      for(ii in seq_len(p-l)){
        INT2 <- INT[CB[,ii], , drop=FALSE] 
        use <- which(colSums(INT2) == Dim-1L)
        for(dd in use){
          inter <- CB[INT2[,dd]==1L,ii] # this will be intersection
          vert <- vert + 1L 
          CCtemp <- cbind(CCtemp, c(inter,k)) 
          lambda <- (U[k]-VTsum_u[ii])/(VTsum_uu[dd]-VTsum_u[ii])
          VTtemp <- cbind(VTtemp,
                          lambda*VT1_u[,dd] + (1-lambda)*VT1_uu[,ii])
        }
      }
    }
  }
  if(any(both)){ 
    for(ll in which(both)){
      vert <- vert + 1L 
      CCtemp <- cbind(CCtemp, CC1[,ll]) 
      VTtemp <- cbind(VTtemp, VT1[,ll]) 
    }
  }
  return(list(VTtemp=VTtemp, CCtemp=CCtemp, vert=vert))
}
