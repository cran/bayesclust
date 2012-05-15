\name{bayesclust-package}
\alias{bayesclust-package}
\alias{bayesclust}
\docType{package}
\title{
 Testing and Searching for Clusters in A Hierarchical Bayes Model
}
\description{
 This package contains a suite of functions that allow 
 the user to carry out the following hypothesis test on genetic data:
 \tabular{l}{
\eqn{H_0} : No clusters \cr 
\eqn{H_1} : 2, 3 or 4 clusters \cr
 }
}
\details{
The hypothesis test is formulated as a model selection problem,
where the aim is to identify the model with the highest posterior probability.  A Hierarchical Bayes
model is assumed for the data.  Note that firstly, the null hypothesis is equivalent to saying that
the population consists of just one cluster.  Secondly, since the functions here only allow the alternative
hypothesis to be either 2, 3 or 4 at any one time, the package allows the user
to test multiple hypotheses while controlling the False Discovery Rate (FDR). 

This is a brief of summary of the test procedure: 
\enumerate{
\item
For a given dataset, compute the empirical posterior probability (EPP) of the null hypothesis
using \code{cluster.test}. EPP will serve as the test statistic in this hypothesis test.
\item
Monitor the convergence of EPP by running \code{plot} on the object returned in Step 1.
\item
Generate the distribution of EPP under the null hypothesis using \code{nulldensity}. This can be 
done concurrent to Steps 1 and 2. Be sure to use the same parameters for Steps 1 and 3 though.
\item
Estimate the \eqn{p}-value of the EPP for this dataset 
using \code{emp2pval}. This function takes the objects returned in 
Steps 1 and 3 as input.
\item
Run \code{cluster.optimal} on significant datasets to pick out optimal clusters.
\item
Run \code{plot} on the object returned in Step 6 to view the optimal clustering/partition of the 
data.
}

For full details on the distributional assumptions, please refer to the papers listed in the references
section.  For further details on the individual functions, please refer to their respective help
pages and the examples.
}
\author{
George Casella \email{casella@stat.ufl.edu} and Claudio Fuentes \email{cfuentes@stat.ufl.edu}
and Vik Gopal \email{viknesh@stat.ufl.edu}

Maintainer: Vik Gopal <viknesh@stat.ufl.edu>
}
\references{
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)

  Gopal, V. and Fuentes, C. and Casella, G. (2012) bayesclust: An R package 
  for Testing and Searching for Significant Clusters.
  \emph{Journal of Statistical Software} \bold{47}(14), 1-21.
  \url{http://www.jstatsoft.org/v47/i14/}
}
\keyword{ package }
\keyword{ htest }
\keyword{ multivariate }
\keyword{ cluster }
\seealso{
 \code{\link{cluster.test}}, \code{\link{cluster.optimal}}, 
 \code{\link{emp2pval}}, \code{\link{nulldensity}} 
}
