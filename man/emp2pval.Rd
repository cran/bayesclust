\name{emp2pval}
\alias{emp2pval}
\alias{emp2pval.cluster.test}
\alias{emp2pval.cluster.test.reps}
\alias{print.emp2pval}
\title{ Convert Empirical Posterior Probability to P-value }
\description{
\code{emp2pval} converts the Empirical Posterior Probability (EPP) computed by
\code{cluster.test} into a frequentist \eqn{p}-value, which can then be used to
assess the significance of the alternative hypothesis.
}
\usage{
\method{emp2pval}{cluster.test}(x, y, ignore = FALSE)
\method{emp2pval}{cluster.test.reps}(x, y, ignore = FALSE)
}
\arguments{
  \item{x}{ \code{x} is an object of class ``cluster.test'' or 
  ``cluster.test.reps'', which is returned when \code{cluster.test} is run on 
  a dataset. }
  \item{y}{ \code{y} is an object of class ``nulldensity'', which is returned when
\code{nulldensity} is run with the necessary parameters provided. }
  \item{ignore}{ \code{ignore} is a logical variable which specifies whether the parameters
in \code{x} and \code{y} should be matched for consistency.  }
}
\details{
  If \code{ignore} is set to FALSE, then the routine will first check to see 
  if the parameters under which the test was run match exactly with the 
  parameters under which the null distribution was generated. If they were, 
  then the EPP's in the ``cluster.test'' object will be converted to a 
  frequentist \eqn{p}-value by checking the ``nulldensity'' object to see 
  which empirical quantile they fall in. If \code{ignore} is set to TRUE, the 
  same EPP to \eqn{p}-value conversion is carried out, but this time without 
  the preliminary check on the parameters.
  
  If replications had been specified when \code{cluster.test} was called, then 
  \code{emp2pval} uses the mean of the estimates of the posterior 
  probabilities.
}
\value{
  \code{emp2pval} returns a list comprising 2 components.
  \item{param }{This component is a copy of the parameters used when running 
    \code{cluster.test} on the dataset. In the case that \code{emp2pval} was 
    run with \code{ignore} set to TRUE, then this could potentially be 
    different than the parameters under which the \code{nulldensity} 
    object was generated.}
  \item{pvals }{A dataframe with the EPP of the dataset and the corresponding 
    frequentist \eqn{p}-values.}
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\author{ Gopal, V. }
\seealso{ \code{\link{cluster.test}} for information on objects of class ``cluster.test''.

\code{\link{nulldensity}} for information on objects of class ``nulldensity''.
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2, nsim=2000)

# Generate corresponding null density object.
null1 <- nulldensity(nsim=100, n=12, p=2, k=2)

# Convert EPP to p-value
emp2pval(test1, null1)

}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
