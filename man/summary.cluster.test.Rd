\name{summary.cluster.test}
\alias{summary.cluster.test}
\alias{summary.cluster.test.reps}
\alias{print.summary.cluster.test}
\alias{print.summary.cluster.test.reps}
\title{ Provides Summary of Cluster Test Result }
\description{
This function provides a summary of objects of class ``cluster.test''.
}
\usage{
\method{summary}{cluster.test}(object, \ldots)
\method{summary}{cluster.test.reps}(object, \ldots)
}
\arguments{
  \item{object}{ \code{object} must be of class ``cluster.test'' or 
  ``cluster.test.reps'', which are the classes of objects returned after 
  running \code{cluster.test}. }
  \item{\dots}{ For later expansion. }
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\author{ Gopal, V. }
\seealso{ 
\code{\link{cluster.test}} for further information on
objects of class ``cluster.test'' and ``cluster.test.reps''.

\code{\link{plot.cluster.test}} and \code{\link{plot.cluster.test.reps}} to 
monitor convergence of computation of posterior probability.}

\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2, nsim=2000, replications=3)

# Summary output
summary(test1)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
