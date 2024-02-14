library(devtools)

devtools::install_github("ohrankwon/eenvlp")
library(eenvlp)

data(cereal)

X = cereal$x
Y = cereal$y

r = dim(Y)[2]
n = dim(Y)[1]

list.prederr = list(fit.eenv = c(),
                    fit.env = c(),
                    fit.ridge = c())
list.fit = list(fit.eenv = c(),
                    fit.env = c(),
                    fit.ridge = c())

for(rep.num in 1:n){
  
  set.seed(rep.num*12+3456789)
  
  # Training/Testing Sets

  X.train <- X[-rep.num, ]
  Y.train <- Y[-rep.num, ]
  X.test <- X[rep.num, ]
  Y.test <- Y[rep.num, ]

  lambda.seq = 10^(seq(-1, -7, length.out=20))
  u.seq = c(0:r)
  ind.cv = sample(n-1, n-1)

  # Choosing u for envelope
  cv.u.env = c()
  for(u.dim in u.seq){
    cv.u.env = c(cv.u.env, eenvlp::cv_eenvlp(X.train, Y.train, u=u.dim, lamb=1e-8, m=n-2, nperm=1, index=ind.cv)) #loocv
  }
  chosen.u.env = u.seq[which.min(cv.u.env)[1]]

  # Choosing u and lamb for enhanced envelope
  cv.lamb.env = matrix(NA,length(lambda.seq),length(u.seq))
  for(i in 1:length(lambda.seq)){
    for(j in 1:length(u.seq)){
      cv.lamb.env[i,j] = eenvlp::cv_eenvlp(X.train, Y.train, u=u.seq[j], lamb=lambda.seq[i], m=n-2, nperm=1, index=ind.cv) #loocv
    }
  }
  chosen = which(cv.lamb.env == min(cv.lamb.env), arr.ind=TRUE)
  chosen.lamb.eenv = lambda.seq[chosen[1]]
  chosen.u.eenv = u.seq[chosen[2]]

  # Estimators
  list.fit[[1]] = eenvlp::eenvlp(X.train, Y.train, u=chosen.u.eenv, lamb=chosen.lamb.eenv) # enhanced envelope
  list.fit[[2]] = eenvlp::eenvlp(X.train, Y.train, u=chosen.u.env, lamb=1e-08) # envelope
  list.fit[[3]] = eenvlp::eenvlp(X.train, Y.train, u=r, lamb=lambda.seq[which.min(cv.lamb.env[,length(u.seq)])]) # ridge
  
  # Prediction
  for(ind.method in c(1:3)){
    Y.pred <- eenvlp::pred_eenvlp(list.fit[[ind.method]], X.test)
    resi <- as.matrix(Y.test - Y.pred)
    #sprederr <- apply(resi, 1, function(x) sum(x ^ 2))
    list.prederr[[ind.method]] <- c(list.prederr[[ind.method]], (sum(resi^2)/r) ) ## mean squared error
  }
  
}

## Prediction errors
print( round(unlist(lapply(list.prederr, mean)),3) )