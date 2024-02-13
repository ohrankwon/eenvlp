eenvlp = function (X, Y, u, lamb, init = NULL, X.scale = TRUE) {
  X <- as.matrix(X)
  X <- scale(X, center = TRUE, scale = X.scale)
  X_center = attributes(X)$`scaled:center`; X_scale = attributes(X)$`scaled:scale`
  Y <- as.matrix(Y)
  
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) 
    stop("u must be an integer between 0 and r.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
    if (nrow(init) != r || ncol(init) != u) 
      stop("The dimension of init is wrong.")
  }
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  svdX <- svd(X)
  betaOLS <- sigYX %*% svdX$v %*% diag(1/(svdX$d^2/n+lamb)) %*% t(svdX$v)
  
  M <- sigY - tcrossprod(betaOLS, sigYX)    
  U <- sigY - M
  tmp <- envMU(M, U, u, initial = init)
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  
  if (u == 0) {
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    muhat <- colMeans(Y)
    betahat <- matrix(0, r, p)
    Sigmahat <- sigY
    #loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
  }
  else if (u == r) {
    etahat <- betaOLS
    Omegahat <- M
    Omega0hat <- NULL
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    betahat <- betaOLS
    Sigmahat <- M
    #loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
  }
  else {
    etahat <- crossprod(Gammahat, betaOLS)
    betahat <- Gammahat %*% etahat
    muhat <- colMeans(Y) - betahat %*% colMeans(X)
    Omegahat <- crossprod(Gammahat, M) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
    Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
    Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, 
                                                  Gamma0hat)
    #loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, 
              beta = betahat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, 
              Omega0 = Omega0hat, #loglik = loglik, 
              n = n, X.center = X_center, X.scale = X_scale))
}



#------------------------------ Other dependent functions

ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
    
  # Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- eigen(X)
  d <- s$values
  nz <- d > tol * d[1]
  structure(
    if(any(nz)) s$vectors[, nz] %*% (t(s$vectors[, nz])/d[nz]) else X,
    dimnames = dnx[2:1])
}
    
GE <- function(A) {
  
  # Gaussian elimination, p must be less than or equal to n
  a <- dim(A)
  n <- a[1]
  p <- a[2]
  idx <- rep(0, p)
  res.idx <- 1:n
  
  i <- 1
  while (i <= p) {
    tmp <- max(abs(A[res.idx, i]))
    Stmp <- setdiff(which(abs(A[, i]) == tmp), idx)
    idx[i] <- Stmp[1]
    res.idx <- setdiff(res.idx, idx[i])
    for (j in 1:(n-i)) {
      A[res.idx[j], ] <- A[res.idx[j], ] - A[res.idx[j], i] / A[idx[i], i] * A[idx[i], ]
    }
    i <- i + 1			
  }
  c(idx, res.idx)
}

envMU <- function(M, U, u, initial = NULL) {
  
  dimM <- dim(M)
  dimU <- dim(U)
  r <- dimM[1]
  
  if (dimM[1] != dimM[2] & dimU[1] != dimU[2]) stop("M and U should be square matrices.")
  if (dimM[1] != dimU[1]) stop("M and U should have the same dimension.")
  #if (qr(M)$rank < r) stop("M should be positive definite.")
  if (u > r & u < 0) stop("u should be between 0 and r.")
  
  
  if (u == 0) {
    
    Gammahat <- NULL
    Gamma0hat <- diag(r)
    MU <- M + U
    tmp.MU <- eigen(MU)
    objfun <- NA #sum(log(tmp.MU$values))
    
  } else if (u == r) {
    
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    tmp.M <- eigen(M)
    objfun <- NA #sum(log(tmp.M$values))
    
  } else if (u == 1) {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])
      
      #		if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2)
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- sum(log(e1$values)) + sum(log(e2$values))
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #			if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- sum(log(e1$values)) + sum(log(e2$values))
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2)
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj4 <- sum(log(e1$values)) + sum(log(e2$values))
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1], ])
    
    i <- 1
    while (i < maxiter) {
      
      fobj <- function(x) {
        T1 <- crossprod(x, x)
        T2 <- crossprod(x, M) %*% x
        T3 <- crossprod(x, invMU) %*% x
        -2 * log(T1) + log(T2) + log(T3)
      }
      
      gobj <- function(x) {
        W1 <- crossprod(x, x)
        T1 <- x / as.vector(W1)
        W2 <- crossprod(x, M) %*% x
        T2 <- M %*% x / as.vector(W2)
        W3 <- crossprod(x, invMU) %*% x
        T3 <- invMU %*% x / as.vector(W3)
        -2 * T1 + T2 + T3
      }
      
      res <- stats::optim(Ginit, fobj, gobj, method = "BFGS")
      g <- as.matrix(res$par)
      a <- qr(g)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))	
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  } else if (u == r - 1 & u != 1) {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
      
      #	  if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- sum(log(e1$values)) + sum(log(e2$values))		
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #	    if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- sum(log(e1$values)) + sum(log(e2$values))			
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2) 
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
      obj4 <- sum(log(e1$values)) + sum(log(e2$values))			
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    #	    }
    #	  }
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])		
    
    j <- GEidx[r]
    
    g <- as.matrix(Ginit[j, ])
    t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
    t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
    
    GUGt2 <- g + t2
    GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
    
    GVGt2 <- g + t3
    GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
    
    invC1 <- chol2inv(chol(GUG))
    invC2 <- chol2inv(chol(GVG))
    
    fobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2	
      T3 <- invC2 %*% tmp3
      -2 * log(1 + sum(x^2)) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
    }
    
    gobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2	
      T3 <- invC2 %*% tmp3
      -4 	* x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	
    }
    
    i <- 1
    while (i < maxiter) {
      
      res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
      Ginit[j, ] <- res$par
      
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))	
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r, drop = FALSE]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  } else {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
      
      #		if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- sum(log(e1$values)) + sum(log(e2$values))		
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #			if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- sum(log(e1$values)) + sum(log(e2$values))			
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2) 
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
      obj4 <- sum(log(e1$values)) + sum(log(e2$values))			
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    #			}
    #		}
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])
    
    
    GUG <- crossprod(Ginit, (M %*% Ginit))	
    GVG <- crossprod(Ginit, (invMU %*% Ginit))		
    
    
    t4 <- crossprod(Ginit[GEidx[(u+1):r],], Ginit[GEidx[(u+1):r], ]) + diag(u)
    i <- 1
    while (i < maxiter) {
      
      for (j in GEidx[(u+1):r]) {
        g <- as.matrix(Ginit[j, ])
        t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
        t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
        
        GUGt2 <- g + t2
        GUG <- GUG - tcrossprod(GUGt2, GUGt2) * M[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
        
        t4 <- t4 - tcrossprod(g, g)
        invC1 <- chol2inv(chol(GUG))
        invC2 <- chol2inv(chol(GVG))
        invt4 <- chol2inv(chol(t4))				
        
        fobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2	
          T3 <- invC2 %*% tmp3
          -2 * log(1 + x %*% T1) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
        }
        
        gobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2	
          T3 <- invC2 %*% tmp3
          -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	
        }
        
        res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
        Ginit[j, ] <- res$par
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        GUGt2 <- g + t2
        GUG <- GUG + tcrossprod(GUGt2, GUGt2) * M[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
        
        
      }
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))	
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  }
  return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat, objfun = objfun))
}
env <- function(X, Y, u, asy = TRUE, init = NULL) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
    if (nrow(init) != r || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigX <- ginv(sigX)
  betaOLS <- sigYX %*% invsigX
  
  U <- tcrossprod(betaOLS, sigYX) 
  M <- sigY - U
  
  tmp <- envMU(M, U, u, initial = init)
  
  
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL 
  
  if (u == 0) {
    
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    muhat <- colMeans(Y)
    betahat <- matrix(0, r, p)
    Sigmahat <- sigY
    loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) ratio <- matrix(1, r, p)
    
  } else if (u == r) {
    
    etahat <- betaOLS
    Omegahat <- M
    Omega0hat <- NULL
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    betahat <- betaOLS
    Sigmahat <- M
    loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- matrix(1, r, p)
    }
    
  } else {
    
    etahat <- crossprod(Gammahat, betaOLS)
    betahat <- Gammahat %*% etahat
    muhat <- colMeans(Y) - betahat %*% colMeans(X)
    Omegahat <- crossprod(Gammahat, M) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
    Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
    Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
      invOmega0hat <- chol2inv(chol(Omega0hat))
      invOmegahat <- chol2inv(chol(Omegahat))
      temp <- kronecker(etahat %*% tcrossprod(sigX, etahat), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
      temp2 <- kronecker(t(etahat), Gamma0hat)
      covMatrix <- kronecker(invsigX, Sigma1) + temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- asyFm / asySE
    }    
  }
  
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, beta = betahat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, ratio = ratio))
}

p.envMU <- function(M, U, u, initial = NULL) {
  
  dimM <- dim(M)
  dimU <- dim(U)
  r <- dimM[1]
  
  if (dimM[1] != dimM[2] & dimU[1] != dimU[2]) stop("M and U should be square matrices.")
  if (dimM[1] != dimU[1]) stop("M and U should have the same dimension.")
  if (qr(M)$rank < r) stop("M should be positive definite.")
  if (u > r & u < 0) stop("u should be between 0 and r.")
  
  
  if (u == 0) {
    
    Gammahat <- NULL
    Gamma0hat <- diag(r)
    MU <- M + U
    tmp.MU <- eigen(MU)
    objfun <- sum(log(tmp.MU$values))
    
  } else if (u == r) {
    
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    tmp.M <- eigen(M)
    objfun <- sum(log(tmp.M$values))
    
  } else if (u == 1) {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])
      
      #		if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- Re( sum(log(eig1$values)) + sum(log(eig2$values)) )
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2)
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- Re( sum(log(e1$values)) + sum(log(e2$values)) )
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #			if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- Re( sum(log(e1$values)) + sum(log(e2$values)) )
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2)
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj4 <- Re( sum(log(e1$values)) + sum(log(e2$values)) )
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1], ])
    
    i <- 1
    while (i < maxiter) {
      
      fobj <- function(x) {
        T1 <- crossprod(x, x)
        T2 <- crossprod(x, M) %*% x
        T3 <- crossprod(x, invMU) %*% x
        -2 * log(T1) + log(T2) + log(T3)
      }
      
      gobj <- function(x) {
        W1 <- crossprod(x, x)
        T1 <- x / as.vector(W1)
        W2 <- crossprod(x, M) %*% x
        T2 <- M %*% x / as.vector(W2)
        W3 <- crossprod(x, invMU) %*% x
        T3 <- invMU %*% x / as.vector(W3)
        -2 * T1 + T2 + T3
      }
      
      res <- stats::optim(Ginit, fobj, gobj, method = "BFGS")
      g <- as.matrix(res$par)
      a <- qr(g)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  } else if (u == r - 1 & u != 1) {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- Re( sum(log(eig1$values)) + sum(log(eig2$values)) )
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
      
      #	  if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- Re( sum(log(eig1$values)) + sum(log(eig2$values)) )
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)	
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #	    if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)		
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2) 
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
      obj4 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)		
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    #	    }
    #	  }
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])		
    
    j <- GEidx[r]
    
    g <- as.matrix(Ginit[j, ])
    t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
    t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
    
    GUGt2 <- g + t2
    GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
    
    GVGt2 <- g + t3
    GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
    
    invC1 <- chol2inv(chol(GUG))
    invC2 <- chol2inv(chol(GVG))
    
    fobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2	
      T3 <- invC2 %*% tmp3
      -2 * log(1 + sum(x^2)) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
    }
    
    gobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2	
      T3 <- invC2 %*% tmp3
      -4 	* x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	
    }
    
    i <- 1
    while (i < maxiter) {
      
      res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
      Ginit[j, ] <- res$par
      
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r, drop = FALSE]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  } else {
    
    maxiter = 100
    ftol = 1e-3
    
    if (!is.null(initial)) {
      init <- initial
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- Re( sum(log(eig1$values)) + sum(log(eig2$values)) )
    } else {
      MU <- M + U
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
      
      midmatrix <- U
      startv <- function(a) t(a) %*% midmatrix %*% a
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
      
      #		if (qr(MU)$rank == r) {
      eig1 <- eigen(t(init) %*% M %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1 <- Re( sum(log(eig1$values)) + sum(log(eig2$values)) )
      
      midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
      tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
      tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
      init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
      e1 <- eigen(t(init.MU) %*% M %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj2 <- Re( sum(log(e1$values)) + sum(log(e2$values)) )		
      if (obj2 < obj1) {
        init <- init.MU
        obj1 <- obj2
      }
      
      #			if (qr(M)$rank == r) {
      tmp.M <- eigen(M)
      midmatrix <- U
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
      obj3 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)		
      if (obj3 < obj1) {
        init <- init.M
        obj1 <- obj3
      }
      
      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
      midmatrix <- invM2 %*% tcrossprod(U, invM2) 
      tmp2.M <- apply(tmp.M$vectors, 2, startv)
      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
      
      e1 <- eigen(t(init.M) %*% M %*% init.M)
      e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
      obj4 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)		
      if (obj4 < obj1) {
        init <- init.M
        obj1 <- obj4
      }
    }
    #			}
    #		}
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])
    
    
    GUG <- crossprod(Ginit, (M %*% Ginit))	
    GVG <- crossprod(Ginit, (invMU %*% Ginit))		
    
    
    t4 <- crossprod(Ginit[GEidx[(u+1):r],], Ginit[GEidx[(u+1):r], ]) + diag(u)
    i <- 1
    while (i < maxiter) {
      
      for (j in GEidx[(u+1):r]) {
        g <- as.matrix(Ginit[j, ])
        t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
        t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
        
        GUGt2 <- g + t2
        GUG <- GUG - tcrossprod(GUGt2, GUGt2) * M[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
        
        t4 <- t4 - tcrossprod(g, g)
        invC1 <- chol2inv(chol(GUG))
        invC2 <- chol2inv(chol(GVG))
        invt4 <- chol2inv(chol(t4))				
        
        fobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2	
          T3 <- invC2 %*% tmp3
          -2 * log(1 + x %*% T1) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
        }
        
        gobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2	
          T3 <- invC2 %*% tmp3
          -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	
        }
        
        res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
        Ginit[j, ] <- res$par
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        GUGt2 <- g + t2
        GUG <- GUG + tcrossprod(GUGt2, GUGt2) * M[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
        
        
      }
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
      obj5 <- Re( sum(log(e1$values)) + sum(log(e2$values))	)
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  }
  return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat, objfun = objfun))
}

