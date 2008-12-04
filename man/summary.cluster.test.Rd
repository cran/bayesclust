\name{summary.cluster.test}
\alias{summary.cluster.test}
\alias{print.summary.cluster.test}
\title{ Provides Summary of Cluster Test Result }
\description{
This function provides a summary of objects of class ``cluster.test''.
}
\usage{
\method{summary}{cluster.test}(object, \ldots)
}
\arguments{
  \item{object}{ \code{object} must be of class ``cluster.test'', which
is the class of object returned after running \code{cluster.test}. }
  \item{\dots}{ For later expansion. }
}
\references{ 
  Fuentes, C. and Casella, G. (2008) "Testing for the Existence of Clusters"
  \url{http://www.stat.ufl.edu/~casella/Papers/paper-v3.pdf}

  Gopal, V. "BayesClust User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/bayesclust/clust.html}
}
\author{ Gopal, V. }
\seealso{ 
\code{\link{cluster.test}} for further information on
objects of class ``cluster.test''.

\code{\link{plot.cluster.test}} to monitor convergence of computation of
posterior probability.}

\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2, nsim=5000)

# Summary output
summary(test1)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
