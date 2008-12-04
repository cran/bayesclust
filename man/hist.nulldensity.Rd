\name{hist.nulldensity}
\alias{hist.nulldensity}
\title{ Plot Histogram of Null Density }
\description{
\code{hist.nulldensity} operates on objects of class 
``nulldensity''. It generates a histogram of the values.
}
\usage{
\method{hist}{nulldensity}(x, \ldots)
}
\arguments{
  \item{x}{ An object of class ``nulldensity''. }
  \item{\dots}{ Additional arguments to be passed to the \code{hist} method. }
}
\references{ 
  Fuentes, C. and Casella, G. (2008) "Testing for the Existence of Clusters"
  \url{http://www.stat.ufl.edu/~casella/Papers/paper-v3.pdf}

  Gopal, V. "BayesClust User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/bayesclust/clust.html}
}
\author{ Gopal, V. }
\seealso{ \code{\link{nulldensity}} for further information on 
objects of class ``nulldensity''.

\code{\link{hist}} for details on arguments that can be passed through \code{\dots}
}
\examples{
# Generate null density object.
null1 <- nulldensity(nsim=100, n=12, p=2, k=2)
hist(null1, main="Null density")

}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
