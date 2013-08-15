\name{coef.glinternet}
\alias{coef.glinternet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Return main effect and interaction coefficients.
}
\description{Returns the actual main effect and interaction coefficients
  that satisfy the sum constraints in a linear interaction model. See
  the paper below for details.
}
\usage{
\method{coef}{glinternet}(object, lambdaIndex = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fitted \code{"glinternet"} model object.}
  \item{lambdaIndex}{Index of \code{lambda} value at which to extract
    coefficients. If NULL, return the coefficients for all values of
    \code{lambda} in the path.}
  \item{\dots}{Not used.}
}
\details{Returns the actual main effect and interaction
  coefficients. These satisfy the sum constraints in the original linear
interaction model.}
\value{A list of length \code{lambda} if \code{lambdaIndex} is not
  provided. Otherwise, of length \code{lambdaIndex}. Each component (for
  each lambda) is
  itself a list, with components
  \item{mainEffects}{A list with components \code{cat} and \code{cont},
    each an index vector of the
     categorical and continuous (resp) main-effect
    variables. Just as in \code{activeSet}, the indexing is separate for
    each type of variable. See ?glinternet for details.}
  \item{mainEffectsCoef}{List of coefficients for the main effects in
    \code{mainEffects}, also with names \code{cat} and \code{cont}}
  \item{interactions}{List of interactions, with components
    \code{contcont}, \code{catcont} and \code{catcat}, each 2-column
    matrices of variable indices.}
  \item{interactionsCoef}{List of interaction coefficients for
    \code{interactions}, also with names \code{contcont}, \code{catcont}
    and \code{catcat}. For categorical-categorical interactions, each
    is provided as a L1 x L2 matrix.}
}
\references{
"Learning interactions via hierarchical group-lasso regularization"
}
\author{Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}}

\seealso{
\code{glinternet.cv}, \code{predict.glinternet},
\code{predict.glinternet.cv}, \code{plot.glinternet.cv}, \code{glinternet}
}

\examples{
Y = rnorm(100)
X = matrix(rnorm(100*10), nrow=100)
numLevels = rep(1, 10)
fit = glinternet(X, Y, numLevels)
coeffs = coef(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}