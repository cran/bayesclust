\name{plot.cluster.test}
\alias{plot.cluster.test}
\title{ Plot Running Posterior Probabilities to Monitor Convergence }
\description{
\code{plot.cluster.test} will plot the running posterior probabilities
that were computed by \code{cluster.test}. These plots can be used to
monitor convergence of the computation carried out by \code{cluster.test}.
}
\usage{
\method{plot}{cluster.test}(x, \ldots)
}
\arguments{
  \item{x}{ \code{x} is an object of class ``cluster.test'', which is
returned when \code{cluster.test} is run on a dataset. }
  \item{\dots}{ Additional arguments to be passed to the inherited \code{plot}
function. }
}
\references{ 
  Fuentes, C. and Casella, G. (2008) "Testing for the Existence of Clusters"
  \url{http://www.stat.ufl.edu/~casella/Papers/paper-v3.pdf}

  Gopal, V. "BayesClust User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/bayesclust/clust.html}
}
\author{ Gopal, V. }
\seealso{ \code{\link{cluster.test}} for further information on
objects of class ``cluster.test''.

\code{\link{plot}} for details on arguments that can be passed through \code{\dots}
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2, nsim=5000)

# Plot the running posterior probabilities to monitor convergence
plot(test1)

}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
