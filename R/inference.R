inference <- function(gfi, v, alpha=0.05){ 
  out <- numeric(4L)
  names(out) <- c("mean","median","low","up")
  vertex <- gfi$VERTEX[v,]
  weight <- gfi$WEIGHT
  out[1] <- sum(vertex*weight) # mean
  h <- cbind(vertex,weight)
  hsort <- h[order(h[,1L]),] # gx.sort(h,1L)
  hsum <- cumsum(hsort[,2L])
  ci_u <- min(which(hsum >= 1-alpha/2)) 
  ci_l <- min(which(hsum >= alpha/2))   
  ci_m <- min(which(hsum >= 0.5))
  out[3] <- hsort[ci_l,1L] # lower bound
  out[4] <- hsort[ci_u,1L] # upper bound
  out[2] <- hsort[ci_m,1L] # estimate (median)
  out
}
