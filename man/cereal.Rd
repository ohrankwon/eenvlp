\name{cereal}
\alias{cereal}

\title{Cereal data}
\usage{
data(cereal)
}

\description{
Cereal data contain 15 cereals on which 145 infrared spectra and two chemical properties, starch and ash, are measured. 
}

\details{The data are from Varmuza and Filzmoser (2016).}

\value{
A list with the following elements:
  \item{x}{A matrix for predictors: 15 rows and 145 columns.}
  \item{y}{A matrix for responses: 15 rows and 2 columns.}
}

\references{
Varmuza, K. and Filzmoser, P. (2016), Introduction to multivariate statistical analysis in chemometrics, CRC press.
}

\examples{
# load data set
data(cereal)

# the dimension of predictors
dim(cereal$x)

# the dimension of responses
dim(cereal$y)
}


\keyword{data set}
