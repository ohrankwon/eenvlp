\name{cv_eenvlp}
\alias{cv_eenvlp}
\title{Cross validation for eenvlp}
\description{
 Compute the prediction error for the enhanced response envelope estimator using m-fold cross validation.
}
\usage{
cv_eenvlp(X, Y, u, lamb, m, nperm, index = NULL, X.scale = TRUE)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{lamb}{A regularization parameter. A positive value.}  
  \item{m}{A positive integer that is used to indicate m-fold cross validation.}
  \item{nperm}{A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.}
  \item{index}{Either NULL or an n vector that is a list of indices from 1 to n. If NULL (the default), then the samples are ordered randomly and split. If a vector is provided, then the samples are ordered accroding to the given list of indices and split.} 
  \item{X.scale}{A boolean. If TRUE (the default), then the columns of training X are divided by their standard deviations. If False, no scaling is done.}
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension u, the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. As Y is multivariate, the identity inner product is used for computing the prediction errors.
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
}

\examples{
data(cereal)
cvPE <- cv_eenvlp(cereal$x, cereal$y, u=1, lamb=1e-07, m=2, nperm=1)
cvPE
}

