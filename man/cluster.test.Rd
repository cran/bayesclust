\name{cluster.test}
\alias{cluster.test}
\title{ Compute Posterior Probabilities for Dataset }
\description{
\code{cluster.test} computes the empirical posterior probability (EPP) of the
null hypothesis in the following test:
\tabular{l}{
\eqn{H_0} : No clusters \cr 
\eqn{H_1} : \eqn{k} clusters \cr
}
where \eqn{k} takes an integer value strictly greater than 1. 
}
\usage{
cluster.test(data, nsim = 500000, aR = 0.4, p = 2, k = 2, 
	a = 2.01, b = 0.990099, tau2 = 1, mcs=0.1, file = "", label = "data",
        replications = 1)
}
\arguments{
  \item{data}{ \code{data} should be a matrix, even when the response is
univariate. The number of rows should equal the number of observations,
and the number of columns should equal to \code{p}. }
  \item{nsim}{ As the Bayes Factor (BF) for the hypothesis is computed
with the aid of a Metropolis-Hastings (MH) MCMC algorithm, the number of
simulations has to be specified. It is recommended that \code{nsim}
be at least 500,000, and that replications be carried out in order to
monitor convergence. }
  \item{aR}{ The candidate distribution in the MH algorithm is
a mixture of \eqn{g} and a random walk. \code{aR}, which must be a value
between 0 and 1, specifies the percentage of time that the random walk is chosen.
Please see the references below for further details on \eqn{g}.}
  \item{p}{ The observations are assumed to come from a multivariate normal
distribution, of length \code{p}. }
  \item{k}{\code{k} must take an integer value strictly greater than 1. It specifies
the alternative hypothesis in the test. }
  \item{a}{ \code{a} is a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{b}{ Like \code{a}, \code{b} is also a hyperparameter for the prior on \eqn{\sigma^2}.
Further details can be found in the references below. }
  \item{tau2}{ \code{tau2} is a hyperparameter for the prior on the mean \eqn{\mu} for
each cluster. }
  \item{mcs}{ \code{mcs} stands for Minimum Cluster Size. It should be a value between 0
and 1. It instructs the test procedure to only consider clusters of a certain minimum
size. }
  \item{file}{ This argument is a character string. If specified, the output object will
be saved to this (binary) file. It can be loaded, inspected and alterered later in
subsequent R sessions using \code{load}. If left unspecified, the object will not
be saved to a file and could be lost on quitting the R session. }
  \item{label}{ \code{label} serves to name the dataset in any given
hypothesis test. }
  \item{replications}{ \code{replications} will instruct the function to 
perform multiple replications of the test, allowing the user to monitor 
convergence of the estimates of the posterior probabilities. }
}
\details{
Since the hypothesis test is carried out in a Bayesian framework, the
Bayes Factor has to be calculated. As this is an integral over a huge space, the sum
is estimated using MCMC. Certain
portions of \code{cluster.test} have been coded in C in order to speed up the simulations.

}
\value{
If \code{replications} is set to 1, then the output from this function is a 
list object consisting of three components. It will
be assigned S3 class ``cluster.test'', and can then serve as input to \code{plot}
or \code{emp2pval}. The components of this list object are:
  \item{param}{ This component exists purely for bookkeeping purposes, in
the sense that these parameters, under which the test was run, will be checked against
the parameters used to generate the null distribution of the test statistic. By default,
conversion of the empirical posterior probability to a \eqn{p}-value will only proceed if the
parameters match. The user does have the option of ignoring this check though. See
\code{emp2pval} for further details.}
  \item{iterations}{ This is a vector of indices that will be used to plot the
running posterior probabilities when \code{plot} is
called. Since it is superfluous to keep the entire chain, the concept of 'thinning'
is applied - only every 500th iteration is stored. Counting begins backwards from
the most recent iteration. }
  \item{ClusterStat}{ This is a 'thinned' vector of running posterior probabilities.
The values that are kept correspond
exactly to those in the preceding \code{iterations} component. When
\code{summary.cluster.test} is called, the final posterior probability is extracted and
printed in a readable format.}
If \code{replications} is set to an integer greater than 1, then the object 
returned is of S3 class ``cluster.test.reps''. This is just a list of 
objects of class ``cluster.test''. This ``cluster.test.reps'' object can also 
be fed directly into \code{emp2pval} and \code{plot}.
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\author{ Gopal, V. }
\seealso{ \code{\link{plot.cluster.test}}, \code{\link{plot.cluster.test.reps}}
 to monitor convergence of computation of posterior probability.

\code{\link{summary.cluster.test}}, \code{\link{summary.cluster.test.reps}} to 
display the computed final posterior probabilities for each dataset run.
}
\examples{
# Generate random 2-variate data
Y <- matrix(rnorm(24), nrow=12)

# Search for optimal partitioning of data into 2 clusters
test1 <- cluster.test(Y, nsim=1000, p=2, replications=3)

# Plot the running posterior probabilities to monitor convergence
plot(test1)

# Generate corresponding null density object.
null1 <- nulldensity(nsim=100, n=12, p=2, k=2)

# Convert EPP to p-value
emp2pval(test1, null1)
}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
