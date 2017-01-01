\name{glinternet.cv}
\alias{glinternet.cv}
\title{Cross-validation for glinternet}
\description{Does k-fold cross validation for glinternet and returns a
  value of \code{lambda}.}
\usage{
glinternet.cv(X, Y, numLevels, nFolds = 10, lambda=NULL, nLambda=50,
lambdaMinRatio=0.01, interactionCandidates=NULL,
screenLimit=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000,
verbose=FALSE, numCores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\code{X} matrix as in \code{glinternet}.}
  \item{Y}{Target \code{Y} as in \code{glinternet}.}
  \item{numLevels}{Number of levels \code{numLevels} as in \code{glinternet}.}
  \item{nFolds}{Number of folds - default is 10.}
  \item{lambda}{\code{lambda} as in \code{glinternet}.}
  \item{nLambda}{\code{nLambda} as in \code{glinternet}.}
  \item{lambdaMinRatio}{\code{lambdaMinRatio} as in \code{glinternet}.}
  \item{interactionCandidates}{\code{interactionCandidates} as in \code{glinternet}.}
  \item{screenLimit}{\code{screenLimit} as in \code{glinternet}.}
  \item{family}{\code{family} as in \code{glinternet}.}
  \item{tol}{\code{tol} as in \code{glinternet}.}
  \item{maxIter}{\code{maxIter} as in \code{glinternet}.}
  \item{verbose}{\code{verbose} as in \code{glinternet}.}
  \item{numCores}{\code{numCores} as in \code{glinternet}.}
}
\details{The \code{lambda} sequence is computed using all the
  data. \code{nFolds} models are fit, each time with one of the folds
  omitted. The error is accumulated, and the average error and standard deviation over the
  folds is computed. The \code{lambda} value that minimizes the average
  error is returned, and a model with this \code{lambda} is fit to the
  full data set.}
\value{
  An object of class \code{glinternet.cv} with the components
  \item{call}{The user function call.}
  \item{glinternetFit}{Glinternet object fitted on the full data using a
    lambda sequence that terminates at \code{lambdaHat}.}
  \item{fitted}{Vector for fitted values (same length as \code{Y}). This
  is from the model fitted at \code{lambdaHat}.}
  \item{activeSet}{\code{activeSet} is a list variables found for the
    model fitted with \code{lambdaHat}.}
  \item{betahat}{Unstandardized coefficients for the variables in
    \code{activeSet}.}
  \item{lambda}{The actual sequence of \code{lambda} values used for the
    cross validation.}
  \item{lambdaHat}{The value of \code{lambda} that minimizes the cv
    error curve.}
  \item{lambdaHat1Std}{The largest value of \code{lambda} that produces
    a cv error that is within 1 standard deviation of the minimum cv
    error. This will always be at least as large as \code{lambdaHat}.}
  \item{cvErr}{The vector of cross validation errors. Same length as
    \code{lambda}.}
  \item{cvErrStd}{Standard deviation for cv errors across the
    \code{nFolds} folds.}
  \item{family}{The response type.}
  \item{numLevels}{Input number of levels for each variable.}
  \item{nFolds}{The number of folds used.}
}

\author{
Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}
}

\seealso{
\code{glinternet}, \code{predict.glinternet},
\code{predict.glinternet.cv}, \code{plot.glinternet.cv}
}
\examples{
Y = rnorm(100)
numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1)
rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
fit = glinternet.cv(X, Y, numLevels, nFolds=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
