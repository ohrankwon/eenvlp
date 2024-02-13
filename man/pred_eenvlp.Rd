\name{pred_eenvlp}
\alias{pred_eenvlp}
\title{Estimation or prediction for rrenv}
\description{
 Perform estimation or prediction under the reduced rank envelope model.
}
\usage{
pred_eenvlp(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from eenvlp.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector or an n by p matrix.}
}
\details{
This function evaluates the reduced rank envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
}

\examples{
data(cereal)
m <- eenvlp(cereal$x[-1,], cereal$y[-1,], u=1, lamb=1e-07)
pred <- pred_eenvlp(m, cereal$x[1,])
pred

}

