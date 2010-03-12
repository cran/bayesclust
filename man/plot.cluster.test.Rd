\name{plot.cluster.test}
\alias{plot.cluster.test}
\alias{plot.cluster.test.reps}
\title{ Plot Running Posterior Probabilities to Monitor Convergence }
\description{
This function will plot the running posterior probabilities
that were computed by \code{cluster.test}. These plots can be used to
monitor convergence of the computation carried out by \code{cluster.test}.
}
\usage{
\method{plot}{cluster.test}(x, \ldots)
\method{plot}{cluster.test.reps}(x, \ldots)
}
\arguments{
  \item{x}{ \code{x} is an object of class ``cluster.test'', or 
  ``cluster.test.reps'', which are returned when \code{cluster.test} is run on 
  a dataset. }
  \item{\dots}{ Additional arguments to be passed to the inherited \code{plot}
function. }
}
\details{
  If only a single replication of the cluster test is carried out, then single 
  plot is displayed on a graphical device. If more than 1 replication had been 
  carried out, then the plots will be arranged in groups of 4, and displayed 
  interactively.
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\author{ Gopal, V. }
\seealso{ \code{\link{cluster.test}} for further information on
objects of class ``cluster.test'' and ``cluster.test.reps''.

\code{\link{plot}} for details on arguments that can be passed through \code{\dots}
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2, nsim=1000, replications=7)

# Plot the running posterior probabilities to monitor convergence
plot(test1)

}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
