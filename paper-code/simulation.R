library(devtools)
devtools::install_github("ohrankwon/eenvlp")
library(eenvlp)
library(MASS)

# Set-ups
args <- commandArgs(trailingOnly = TRUE)
n = as.numeric(args[1])
p = as.numeric(args[2])
rho = as.numeric(args[3])

rate = p/n
num.iter = 100

# True beta, SigmaX, gamma, eta, .... 
D = diag(c(10,8,2))
alpha = pi/4
beta = pi/4
gamma = pi/4
v1 = c(cos(alpha)*cos(gamma)-sin(alpha)*sin(beta)*sin(gamma),cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma),cos(beta)*sin(gamma))
v2 = c(-sin(alpha)*cos(beta),cos(alpha)*cos(beta),-sin(beta))
v3 = c(-cos(alpha)*sin(gamma)-sin(alpha)*sin(beta)*cos(gamma),cos(alpha)*sin(beta)*cos(gamma)-sin(alpha)*sin(gamma),cos(beta)*cos(gamma))
V = matrix(c(v1,v2,v3),3,3)
Vin = solve(V)
Sigma=V%*%D%*%Vin

e.vec=eigen(Sigma)$vectors
gamma = e.vec[,2:3]
set.seed(11111); eta = matrix(rnorm(2*p),2,p)
eta = t(t(eta)/sqrt(sum(eta^2))*sqrt(10))
beta = gamma%*%eta
SigmaX = rho^abs(outer(1:p,1:p,"-"))

# Generating X
set.seed(19910206)
if(rho==0){
  X = matrix(rnorm(n*p),n,p)
} else {
  X = mvrnorm(n, rep(0,p), SigmaX)
}

# Function to calculate prediction risk and bias -- you need to know true beta and SigmaX
risk_eenvlp <- function(m, true.beta, true.SigmaX) {

  Xmean = m$X.center
  if(is.null(m$X.scale)) {Xsd = rep(1,length(Xmean))} else {Xsd = m$X.scale}
   
  r <- ncol(m$Sigma)
  n <- m$n
  if (is.null(m$Gamma)) {
    u <- 0
  } else {
    u <- ncol(m$Gamma)
  }
   
  if (u == 0) {
    alp.true <- t( - true.beta)
  } else {
    alp.true <- diag(1/Xsd) %*% t(m$beta) - t(true.beta)
  }
   
  alpha.true <- t(alp.true) %*% true.SigmaX %*% alp.true
 
  output = list()
  output[[1]] = sum(diag(alpha.true))
  output[[2]] = t(alp.true)
 
  return( output )
}

# Results (risks, chosen u, bias) will be saved 
result = list()
result[[1]] = matrix(NA,num.iter,7)
result[[2]] = list(n=n, p=p, p_n=p/n, rho=rho)

cum_bias = list()
for(i in 1:5){
  cum_bias[[i]] = matrix(0,3,p)
}

#--------------------

for(seed in 1:num.iter){

  # Generating Y given X
  set.seed(seed)
  e=mvrnorm(n, c(0,0,0), Sigma)
  Y=X%*%t(beta)+e
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  r <- dim(Y)[2]
  
  # Finding tuning parameters using 10-fold cv 
  cv.index = sample(n,n)
  lambda.seq = 10^(seq(1, -1, length.out=100))

  cv.pe.lamb.0 = c() # list of pe when u=0
  for(i in 1:length(lambda.seq)){
    cv.pe.lamb.0 = c(cv.pe.lamb.0, eenvlp::cv_eenvlp(X, Y, u=0, lamb=lambda.seq[i], m=10, nperm=1, index=cv.index))
  }  
  
  cv.pe.lamb.1 = c() # list of pe when u=1
  for(i in 1:length(lambda.seq)){
    cv.pe.lamb.1 = c(cv.pe.lamb.1, eenvlp::cv_eenvlp(X, Y, u=1, lamb=lambda.seq[i], m=10, nperm=1, index=cv.index))
  }  
  
  cv.pe.lamb.2 = c() # list of pe when u=2
  for(i in 1:length(lambda.seq)){
    cv.pe.lamb.2 = c(cv.pe.lamb.2, eenvlp::cv_eenvlp(X, Y, u=2, lamb=lambda.seq[i], m=10, nperm=1, index=cv.index))
  }  
  
  cv.pe.lamb.3 = c() # list of pe when u=3
  for(i in 1:length(lambda.seq)){
    cv.pe.lamb.3 = c(cv.pe.lamb.3, eenvlp::cv_eenvlp(X, Y, u=r, lamb=lambda.seq[i], m=10, nperm=1, index=cv.index))
  }  

  # -- u and lamb chosen for enhanced envelope
  cv.pe.lamb.comb = rbind(cv.pe.lamb.0, cv.pe.lamb.1, cv.pe.lamb.2, cv.pe.lamb.3)
  cv.pe.eenv.min = which(cv.pe.lamb.comb==min(cv.pe.lamb.comb), arr.ind=T)

  # -- u chosen for envelope
  cv.pe.org = c()
  if(rate >= 0.9){small.lamb=1e-08}else{small.lamb=0} # each fold becomes n<p
  for(i in 0:r){
    cv.pe.org = c(cv.pe.org, eenvlp::cv_eenvlp(X, Y, u=i, lamb=small.lamb, m=10, nperm=1, index=cv.index))
  }
  cv.pe.env.min = which(cv.pe.org==min(cv.pe.org))

  # Prediction risk calculation 
  if(rate >= 1.0){small.lamb=1e-08}else{small.lamb=0}  
  PE_Eenv_u_fix = risk_eenvlp(eenvlp::eenvlp(X,Y,u=2,lamb=lambda.seq[which.min(cv.pe.lamb.2)]), true.beta=beta, true.SigmaX=SigmaX)
  PE_Env_u_fix = risk_eenvlp(eenvlp::eenvlp(X,Y,u=2,lamb=small.lamb), true.beta=beta, true.SigmaX=SigmaX)
  PE_Eenv_u_chosen = risk_eenvlp(eenvlp::eenvlp(X,Y,u=cv.pe.eenv.min[1]-1,lamb=lambda.seq[cv.pe.eenv.min[2]]), true.beta=beta, true.SigmaX=SigmaX)
  PE_Env_u_chosen = risk_eenvlp(eenvlp::eenvlp(X,Y,u=cv.pe.env.min[1]-1,lamb=small.lamb), true.beta=beta, true.SigmaX=SigmaX)
  PE_Ridge = risk_eenvlp(eenvlp::eenvlp(X,Y,u=r,lamb=lambda.seq[which.min(cv.pe.lamb.3)]), true.beta=beta, true.SigmaX=SigmaX)

  result[[1]][seed,] = c(
    PE_Eenv_u_fix[[1]],
    PE_Env_u_fix[[1]],
    PE_Eenv_u_chosen[[1]],
    PE_Env_u_chosen[[1]],
    PE_Ridge[[1]],
    cv.pe.eenv.min[1]-1,
    cv.pe.env.min[1]-1
  )

  cum_bias[[1]] =  cum_bias[[1]] + PE_Eenv_u_fix[[2]]
  cum_bias[[2]] =  cum_bias[[2]] + PE_Env_u_fix[[2]]
  cum_bias[[3]] =  cum_bias[[3]] + PE_Eenv_u_chosen[[2]]
  cum_bias[[4]] =  cum_bias[[4]] + PE_Env_u_chosen[[2]]
  cum_bias[[5]] =  cum_bias[[5]] + PE_Ridge[[2]]

}

# Combining results
colnames(result[[1]]) = c("risk_Eenv-u_fix", "risk_Env-u_fix", "risk_Eenv-u_chosen", "risk_Env-u_chosen", "risk_Ridge", "u_Eenv-u_chosen", "u_Env-u_chosen")

meanresult1 = colMeans(result[[1]])
meanresult2 = round(meanresult1,2)
seresult1 = apply(result[[1]], 2, sd)/sqrt(num.iter)
seresult2 = round(seresult1,2)
result[[1]] = rbind(result[[1]],meanresult1,seresult1,ave=meanresult2,se=seresult2)

svdSigmaX <- svd(SigmaX)
sqrt.true.SigmaX <- svdSigmaX$u %*% diag(sqrt(svdSigmaX$d)) %*% t(svdSigmaX$v)

for(i in 1:5){
  bias2 = sqrt.true.SigmaX %*% t(cum_bias[[i]]/num.iter)
  cum_bias[[i]] = sum(bias2^2)
}
bias2 = unlist(cum_bias)
names(bias2) = c("Bias2_Eenv-u_fix", "Bias2_Env-u_fix", "Bias2_Eenv-u_chosen", "Bias2_Env-u_chosen", "Bias2_Ridge")


# Printing results
cat("\n General information \n")
print(unlist(result[[2]]))

cat("\n Prediction risk \n")
print(tail(result[[1]][,1:5],2))

cat("\n Freq table of chosen u from enhanced env")
print(table(result[[1]][1:num.iter,6])) 

cat("\n Freq table of chosen u from env")
print(table(result[[1]][1:num.iter,7])) 

cat("\n Bias Squared \n")
print(round(unlist(bias2),2))
