\name{nulldensity}
\alias{nulldensity}
\title{ Generate Null Distribution of Empirical Posterior Probability }
\description{
In testing the following the hypothesis,
\tabular{l}{
\eqn{H_0} : No clusters \cr
\eqn{H_1} : \eqn{k} clusters \cr
}
\code{nulldensity} generates random variables from the distribution of
the Empirical Posterior Probability (EPP) under the null hypothesis.
}
\usage{
nulldensity(nsim, n, k, mcs=0.2, a=2.01, b=0.990099,
	tau2=1, prop=0.25, p, file="")
}
\arguments{
  \item{nsim}{ This denotes the number of random variables to generate
from the distriution of EPP under the null. It is recommended to
be at least 100,000 when the intention is to carry out multiple testing.
Otherwise 8,000-10,000 iterations will suffice.}
  \item{n}{ \code{n} is the number of observations in the dataset for
testing this hypothesis. See \code{\link{cluster.test}} for more details. }
  \item{k}{ \code{k} specifies the alternative hypothesis being tested. It
must take an integer value strictly greater than 1. }
  \item{mcs}{ \code{mcs} stands for Minimum Cluster Size. It should be a value between 0
and 1. It instructs the test procedure to only consider clusters of a certain minimum
size. }
  \item{a}{ \code{a} is a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{b}{ Like \code{a}, \code{b} is also a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{tau2}{ \code{tau2} is a hyperparameter for the prior on the mean \eqn{\mu} for
each cluster. }
  \item{prop}{ \code{prop} specifies what fraction of the space of partitions under the
null hypothesis should be sampled. It is recommended to be at least 0.25. }
  \item{p}{ The observations are assumed to come from a multivariate normal
distribution, of length \code{p}. }
  \item{file}{ This argument is a character string. If specified, the output object will
be saved to this (binary) file. It can be loaded, inspected and alterered later in
subsequent R sessions using \code{load}. If left unspecified, the object will not
be saved to a file and could be lost on quitting the R session. }
}
\details{
The test statistic (EPP) is computed by the function \code{cluster.test}. In order
to assess the significance of the statistic, it is necessary to obtain the
frequentist \eqn{p}-value of
the calculated statistic. This package achieves this task by simulating the null distribution
of the test statistic with \code{nulldensity} and then extracting the sample quantile using
\code{emp2pval}.

A very small portion of the code has been written in C. The code becomes slower as
\code{k} gets larger in the alternative hypothesis.

For a particular dataset, this function can be run in parallel with \code{cluster.test}.
}
\value{
The object returned is of class ``nulldensity''. It is a list comprising two components.
  \item{param }{ This component, again, exists purely for bookkeeping purposes.
When \code{emp2pval} is called, it takes two mandatory arguments - one of class ``cluster.test''
and the other of class ``nulldensity''. Both these objects have a parameter component, which
should match for the \eqn{p}-value conversion to proceed. }
  \item{gen.values }{ This is a vector of length \code{nsim}, consisting of the simulations
from the null distribution. }
}
\references{ 
  Fuentes, C. and Casella, G. (2008) "Testing for the Existence of Clusters"
  \url{http://www.stat.ufl.edu/~casella/Papers/paper-v3.pdf}

  Gopal, V. "BayesClust User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/bayesclust/clust.html}
}
\author{ Fuentes, C. and Gopal, V. }
\seealso{
\code{\link{cluster.test}} for further information on
objects of class ``cluster.test''.

\code{\link{hist.nulldensity}} which allows the user to plot a histogram of
simulated values in order to view the shape of the null distribution.
}
\examples{
# Generate null density object.
null1 <- nulldensity(nsim=100, n=12, p=2, k=2)
hist(null1)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
