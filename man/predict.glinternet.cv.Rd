\name{predict.glinternet.cv}
\alias{predict.glinternet.cv}
\title{Make predictions from a "glinternetCV" object.}
\description{
Similar to other predict methods, this function returns fitted values on
the response scale. Also gives the option to return the link function.
}
\usage{
\method{predict}{glinternet.cv}(object, X, type = c("response", "link"),
lambdaType=c("lambdaHat", "lambdaHat1Std"), ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{"glinternetCV"} object.}
  \item{X}{Matrix of new values for which to make predictions. Must have
  the same number of variables as during training the model, and all the
  variables must have the same number of levels.}
  \item{type}{Return fitted \code{Y}
    values or the link function.}
  \item{lambdaType}{\code{lambdaHat} corresponds to the lambda value
    that gives the minimum CV error. \code{lambdaHat1Std} picks the
    largest value of lambda for which the CV error is within 1 standard
    error of the minimum.}
  \item{\dots}{Not used.}
}
\details{
Makes predictions using the model fitted at the appropriate lambda value.
}
\value{A vector of predicted values.}

\author{
Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{glinternet}, \code{glinternet.cv}, \code{predict.glinternet}
}
\examples{
Y = rnorm(100)
numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1)
rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
fit = glinternet.cv(X, Y, numLevels, nFolds=3)
max(abs(fit$fitted - predict(fit, X)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
