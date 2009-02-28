\name{fdr.test}
\alias{fdr.test}
\title{ For Testing Multiple Hypotheses, Controlling for FDR }
\description{
Allows the experimenter to assess the significance of clusters in multiple
datasets whilst controlling for the False Discovery Rate (FDR).
}
\usage{
fdr.test(namelist, q=0.05)
}
\arguments{
  \item{namelist}{ The initial argument should be a character vector of names
of objects of class ``emp2pval''. Objects of this class are returned by the function
\code{emp2pval}. }
  \item{q}{ \code{q} should be a value between 0 and 1. It is the maximum FDR
that the experimenter is willing to tolerate.}
}
\details{
  This function implements the FDR controlling procedure
described in Benjamini and Hochberg (1995). This routine will look for the \emph{largest} \eqn{i} such
that \eqn{P_{(i)} \leq (i/m) * q }{P_{(i)} <= (i/m) * q }.
Please refer to the original paper for the notation and details.
}
\value{
  The output simply lists the datasets for which significant clusters were found.
}
\references{
  Benjamini, Y. and Hochberg, Y. (1995)
  \emph{Controlling the False Discovery Rate: A Practical and Powerful Approach to
  Multiple Testing}.
  JRRS (B) Vol. 57, No. 1
}
\author{ Gopal, V. }
\seealso{
\code{\link{emp2pval}} for details on objects of class \code{emp2pval}.
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)
Z <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, p=2)
test2 <- cluster.test(Z, p=2)

# Generate corresponding null density object.
null1 <- nulldensity(nsim=100, n=12, p=2, k=2)

# Convert EPP to p-value
test1.pval <- emp2pval(test1, null1)
test2.pval <- emp2pval(test2, null1)

# Test for significance, controlling for FDR
fdr.test(c("test1.pval", "test2.pval"), q=0.3)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
