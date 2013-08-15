\name{predict.glinternet}
\alias{predict.glinternet}
\title{Make predictions from a "glinternet" object.}
\description{
Similar to other predict methods, this function returns fitted values on
the response scale. Also gives the option to return the link function.
}
\usage{
\method{predict}{glinternet}(object, X, type = c("response", "link"), lambda=NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fitted \code{"glinternet"} model object.}
  \item{X}{Matrix of new values for which to make predictions. Must have
  the same number of variables as during training the model, and all the
  variables must have the same number of levels.}
  \item{type}{Return fitted \code{Y}
    values or the link function.}
  \item{lambda}{User input lambda sequence. Must be subset of the
    lambdas used in fitting. If NULL (the default), predict at all the
    lambdas used during fitting.}
  \item{\dots}{Not used. Other arguments to predict.}
}
\details{
If lambda is not specified, makes predictions at all the fitted lambda values. Users may provide their own lambda sequence, but this must be a subset of the values used to fit the models. 
}
\value{A matrix of predicted values, with columns corresponding to each
  fitted model.}

\author{
Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{glinternet}, \code{glinternet.cv}, \code{predict.glinternet.cv}, \code{plot.glinternet.cv}
}
\examples{
Y = rnorm(100)
numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1)
rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
fit = glinternet(X, Y, numLevels)
max(abs(fit$fitted - predict(fit, X)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
