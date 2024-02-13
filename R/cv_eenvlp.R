cv_eenvlp = function(X, Y, u, lamb, m, nperm, index=NULL, X.scale = TRUE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  prederr <- rep(0, nperm)
  
  for (i in 1:nperm)	{
    if(is.null(index)){ id <- sample(n, n) } else { id <- index }
    Xn <- X[id, ]
    Yn <- Y[id, ]
    Yn <- as.matrix(Yn)
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n / m) + 1) : ceiling(j * n / m)
      id.train <- setdiff(1:n, id.test)
      
      X.train <- Xn[id.train, ]; X.train = scale(X.train, center=TRUE, scale=X.scale)
      Y.train <- Yn[id.train, ]
      X.test <- Xn[id.test, ]
      if (X.scale) {
        X.test = scale(X.test, center=attributes(X.train)$`scaled:center`, scale=attributes(X.train)$`scaled:scale`)
      } else {
        X.test = scale(X.test, center=attributes(X.train)$`scaled:center`)
      }
      Y.test <- Yn[id.test, ]
      n.test <- length(id.test)
      fit <- eenvlp(X.train, Y.train, u, lamb, X.scale = FALSE)
      betahat <- t(fit$beta)
      muhat <- fit$mu
      resi <- as.matrix(Y.test - matrix(1, n.test, 1) %*% t(muhat) - as.matrix(X.test) %*% betahat)
      sprederr <- apply(resi, 1, function(x) sum(x ^ 2))
      prederr[i] <- prederr[i] + sum(sprederr)
    }
  }
  return(sqrt(mean(prederr / n)))
}