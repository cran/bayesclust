\name{combine}
\Rdversion{1.1}
\alias{combine}
\title{
  Combine posterior probabilities of null hypothesis under several tests
}
\description{
  In the case where the user knows the maximum number of clusters, this 
  function allows for the tester to combine several simple hypothesis tests 
  on the same data set to gauge the posterior probabilities of k=2, 3, 4 
  clusters, given the observed data.
}
\usage{
combine(...)
}
\arguments{
  \item{\dots}{
  The arguments to \code{combine} should be of class ``cluster.test'' or 
  ``cluster.test.reps''. These should be the results of running a single 
  dataset through the \code{cluster.test} function with the exact same 
  parameters, but with different settings for \code{k}, i.e. this function 
  is appropriate when testing for differing number of clusters on the same 
  data set.
}
}
\details{
  The objects fed into this function need not all be of the same class. Some 
  could be of class ``cluster.test'' while others could be of class 
  ``cluster.test.reps''. The function will take the mean of the replications 
  when computing the Bayes Factor for each test. For details on the  
  calculations, please refer to the references below.
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
testK2 <- cluster.test(Y, k=2, p=2, nsim=2000, replications=3)
testK3 <- cluster.test(Y, k=3, p=2, nsim=2000)
testK4 <- cluster.test(Y, k=4, p=2, nsim=2000)

combine(testK2, testK3, testK4)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
