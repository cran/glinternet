\name{plot.glinternet.cv}
\alias{plot.glinternet.cv}
\title{Plot CV error from \code{glinternetCV} object.}
\description{Plots the cross validation error against the lambda index. Uses \code{ggplot2} if found on the user's system.
}
\usage{
\method{plot}{glinternet.cv}(x, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{"glinternetCV"} object.}
  \item{\dots}{Not used.}
}

\value{A plot of CV error.}

\author{
Michael Lim and Trevor Hastie\cr Maintainer: Michael Lim \email{michael626@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{glinternet}, \code{glinternet.cv}, \code{predict.glinternet.cv}, \code{predict.glinternet}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ interactions }% __ONLY ONE__ keyword per line
\keyword{group-lasso}
