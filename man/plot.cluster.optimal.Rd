\name{plot.cluster.optimal}
\alias{plot.cluster.optimal}
\title{ Plot the Optimal Clusters Found in Datasets }
\description{
 \code{plot.cluster.optimal} allows the user to plot the
 optimal clusters found by the search algorithm in
 \code{cluster.optimal}.
}
\usage{
\method{plot}{cluster.optimal}(x, varToPlot=c(1,2), clustToDisp=c(1,2,3,4), 
  \ldots)
}
\arguments{
  \item{x}{ \code{x} must be an object of class ``cluster.optimal'', which is
returned when \code{cluster.optimal} is run on a dataset. }
  \item{varToPlot}{ \code{varToPlot} must be vector of length 2, specifying 
which of the \code{p} variables should be plotted. }
  \item{clustToDisp}{ \code{clustToDisp} must be a vector of length at most 4, 
indicating how many of the clusters to display. }
  \item{\dots}{ Additional arguments to be passed to the inherited \code{plot} }
}
\details{
 This function will plot the optimal clusters found, using different colors
to represent the different clusters. The function will plot two variables at a 
time; these can be specified via the \code{varToPlot} argument. If only the 
top cluster had been kept, then this will be plotted. If more than one 
cluster had been retained, then the user can specify which of them should be 
displayed through the \code{clustToDisp} argument. 

 However, for the purpose of keeping the plots readable, only 4 clusterings 
will be plotted on a given window. For the remaining clusterings, the user can 
call the function again with a different argument vector for 
\code{clustToDisp}. 

}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
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
search1 <- cluster.optimal(Y, p=2, keep=5, nsim=1000)

# Plot the best cluster found during search
plot(search1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
