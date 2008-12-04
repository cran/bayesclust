\name{plot.cluster.optimal}
\alias{plot.cluster.optimal}
\title{ Plot the Optimal Clusters Found in Datasets }
\description{
 \code{plot.cluster.optimal} allows the user to plot the
 optimal clusters found by the search algorithm in
 \code{cluster.optimal}.
}
\usage{
\method{plot}{cluster.optimal}(x, \ldots)
}
\arguments{
  \item{x}{ \code{x} must be an object of class ``cluster.optimal'', which is
returned when \code{cluster.optimal} is run on a dataset. }
  \item{\dots}{ Additional arguments to be passed to the inherited \code{plot} }
}
\details{
 This function will plot the optimal cluster found, using different colors
to represent the different clusters. If the experimenter wishes to plot
the 2nd best, or 3rd best cluster, he will have to do this manually.

  This function only works for the cases when the observations in the data are
vectors of length \eqn{p}, where \eqn{p} takes one of the values 2, 3, 4, 5 or 6.
In each of these cases pairwise plots will be generated, one for each dataset.

In the event that the experimenter only wishes to plot certain variables against another,
he/she will have to do this manually, and may wish to refer to the examples below.
}
\references{ 
  Fuentes, C. and Casella, G. (2008) "Testing for the Existence of Clusters"
  \url{http://www.stat.ufl.edu/~casella/Papers/paper-v3.pdf}

  Gopal, V. "BayesClust User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/bayesclust/clust.html}
}
\author{ Gopal, V. }
\seealso{ \code{\link{cluster.optimal}} for further information on
objects of class ``cluster.optimal''.

\code{\link{plot}} for details on arguments that can be passed through \code{\dots}
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
search1 <- cluster.optimal(Y, p=2, keep=5)

# Plot the best cluster found during search
plot(search1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
