#' @importFrom stats runif
fid_sample <- function(VT2, VTsum, U, L){
  high <- VT2>0
  low <- VT2<0
  zero <- VT2==0
  if((any(high) && any(low)) || any(zero)){
    UU <- sign(U-VTsum)
    LL <- sign(L-VTsum)
    SS <- sign(VT2)
    zeronot <- !zero
    if(sum(zero)==length(VT2)){
      MAX <- Inf
      MIN <- -Inf
      temp <- any(UU>0)*any(LL<0)
    }else{
      c1 <- !any(UU[zero]==-1) && !any(LL[zero]==-1)
      c2 <- !any(UU[zero]==1) && !any(LL[zero]==1)
      d1 <- !any(SS==-1)
      d2 <- !any(SS==1)
      if((d1 && c1) || (d2 && c2)){
        Us <- (U-VTsum[zeronot])/VT2[zeronot]
        Ls <- (L-VTsum[zeronot])/VT2[zeronot]
        MAX <- Inf
        MIN <- min(Us,Ls)
        temp <- 1-(atan(MIN)/pi+0.5)
      }else if((d2 && c1) || (d1 && c2)){
        Us <- (U-VTsum[zeronot])/VT2[zeronot]
        Ls <- (L-VTsum[zeronot])/VT2[zeronot]
        MAX <- max(Us,Ls)
        MIN <- -Inf
        temp <- atan(MAX)/pi+0.5 				
      }else{
        HUs <- (U-VTsum[high])/VT2[high]
        HLs <- (L-VTsum[high])/VT2[high]
        Hmax <- max(HUs,HLs)
        Hmin <- min(HUs,HLs)
        LUs <- (U-VTsum[low])/VT2[low]
        LLs <- (L-VTsum[low])/VT2[low]
        Lmax <- max(LUs,LLs)
        Lmin <- min(LUs,LLs)
        if(round(Lmin-Hmax,12)>=0){
          bpos <- -Inf
          tpos <- Hmax
          bneg <- Lmin
          tneg <- Inf
        }else if(round(Hmin-Lmax,12)>=0){
          bpos <- Hmin
          tpos <- Inf
          bneg <- -Inf
          tneg <- Lmax
        }else{
          bpos <- -Inf
          tpos <- Inf
          bneg <- -Inf
          tneg <- Inf
        }
        if(tpos==Inf){
          Pprob <- 1-(atan(bpos)/pi+0.5)
        }else{
          Pprob <- atan(tpos)/pi+0.5
        }
        if(tneg==Inf){
          Nprob <- 1-(atan(bneg)/pi+0.5)
        }else{
          Nprob <- atan(tneg)/pi+0.5
        }
        temp <- Pprob+Nprob
        Pprob <- Pprob/temp
        Nprob <- 1-Pprob
        if(runif(1L) <= Nprob){
          MIN <- bneg
          MAX <- tneg
        }else{
          MIN <- bpos
          MAX <- tpos
        }
      }
    }
    y <- atan(MAX)/pi+0.5
    x <- atan(MIN)/pi+0.5
    u <- x+(y-x)*runif(1L)
    ZZ <- tan(pi*(u-0.5))
    ZZ2 <- ZZ*ZZ
    wt <- exp(-ZZ2/2)*(1+ZZ2)*temp
  }else{
    Us <- (U-VTsum)/VT2
    Ls <- (L-VTsum)/VT2
    MAX <- max(Us,Ls)
    MIN <- min(Us,Ls)
    y <- atan(MAX)/pi+.5; #cdf
    x <- atan(MIN)/pi+.5
    u <- x+(y-x)*runif(1L)
    ZZ <- tan(pi*(u-.5)) #Inverse cdf
    ZZ2 <- ZZ*ZZ
    wt <- exp(-ZZ2/2)*(1+ZZ2)*(y-x)
  }
  return(c(ZZ=ZZ,wt=wt))
}
