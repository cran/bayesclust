\name{cluster.optimal}
\alias{cluster.optimal}
\title{ Search for Optimal Clustering of Dataset }
\description{
  \code{cluster.optimal} will search for the optimal
\eqn{k}-clustering of the dataset.
}
\usage{
cluster.optimal(data, nsim = 100000, aR = 0.4, p = 2, k = 2,
	a = 2.01, b = 0.990099, tau2 = 1, mcs = 0.1, 
	file = "", label = "data", keep = 4)
}
\arguments{
  \item{data}{ \code{data} should be a matrix, even when the response is
univariate. The number of rows should equal the number of observations,
and the number of columns should equal to \code{p}. }
  \item{nsim}{ The algorithm is based on a stochastic search through the
partition space, and \code{nsim} corresponds to the number of points in the
partition space to inspect. It is recommended that \code{nsim} be at least
500,000.}
  \item{aR}{ The Metropolis search algorithm samples points from the space of
partitions according to a mixture of \eqn{g} and a
random walk. \code{aR}, which must be a value between 0 and 1, specifies the
percentage of time that the random walk is chosen. Please see the references
below for further details on the density \code{g}.  }
  \item{p}{ The observations are assumed to come from a multivariate normal
distribution, of length \code{p}. }
  \item{k}{ \code{k} must take an integer value strictly greater than 1. It instructs the algorithm
to search for the optimal partition of the data into \code{k} clusters.}
  \item{a}{ \code{a} is a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{b}{ Like \code{a}, \code{b} is also a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{tau2}{ \code{tau2} is a hyperparameter for the prior on the mean \eqn{\mu} for
each cluster. }
  \item{keep}{ This argument instructs the algorithm to store the top \code{keep} number of
clusters that it finds during it's run. By default, the best 4 clusters that are found will
be kept. }
  \item{mcs}{ \code{mcs} stands for Minimum Cluster Size. It should be a value between 0
and 1. It instructs the algorithm to only consider clusters of a certain minimum
size. }
  \item{file}{ This argument is a character string. If specified, the output object will
be saved to this (binary) file. It can be loaded, inspected and altered later in
subsequent R sessions using \code{load}. If left unspecified, the object will not
be saved to a file and could be lost on quitting the R session. }
  \item{label}{ \code{label} serves to name the dataset in any given
hypothesis test. }
}
\details{
  A Metropolis search algorithm is run to maximise the marginal of \eqn{Y}, that is, \eqn{m(Y | \omega_k)}
where \eqn{\omega_k} is a particular partitioning of the data into
\eqn{k} clusters. The partitions that yield the highest marginal will be deemed to be optimal.
 
  Since the sample space is so large, the algorithm is started at an intelligent
starting point, by running kmeans, but with a crude imposition of the minimum cluster size.
}
\value{
  The object returned is a list consisting of 2 components.
  \item{param}{ The purpose of this component is to store the parameters under
which the algorithm was run.  }
  \item{data }{This component will
contain a table of the best \code{keep} clusters and the computed value of the
marginal corresponding to those clusters. The latter values are kept as a means of assessing
the relative merit of the clusters in the table. }
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\author{ Gopal, V. }
\seealso{ \code{\link{plot.cluster.optimal}} to plot the clustered data.
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
search1 <- cluster.optimal(Y, nsim=1000, p=2, k=4, keep=5)

# Plot the best cluster found during search
plot(search1)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
