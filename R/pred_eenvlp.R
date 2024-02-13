pred_eenvlp <- function(m, Xnew) {
  

  r <- ncol(m$Sigma)
  n <- m$n
  if (is.null(m$Gamma)) {
    u <- 0
  } else {
    u <- ncol(m$Gamma)
  }
  
  if( is.vector(Xnew) ) {
    Xnew <- matrix(Xnew, 1, length(Xnew))
  }
  Xnew <- as.matrix(Xnew)
  
  if( is.null(m$X.scale) ) {
    Xnew <- scale(Xnew, center = m$X.center)
  } else {
    Xnew <- scale(Xnew, center = m$X.center, scale = m$X.scale)
  }
  
  if (u == 0) {
    value <- matrix(1, dim(Xnew)[1], 1) %*% t(m$mu)
    
  } else {
    value <- matrix(1, dim(Xnew)[1], 1) %*% t(m$mu) + as.matrix(Xnew) %*% t(m$beta)
  }
  
  return(
    value
  )
  
}
