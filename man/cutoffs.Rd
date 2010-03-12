\name{cutoffs}
\alias{cutoffs}
\docType{data}
\title{ Table of Cut-Off Points }
\description{
  This is a table of pre-computed cut-off points for testing
  significance of clusters at the \eqn{alpha}=0.05 and 0.01 level.
}
\usage{data(cutoffs)}
\format{
  A table of cut-off points obtained by generating the null distribution
  of the posterior probability using the following parameters:
  \describe{
    \item{\code{n}}{The number of observations in the dataset.}
    \item{\code{mcs}}{\code{mcs} stands for Minimum Cluster Size. }
    \item{\code{p}}{The length of the vector of each observation. For example, \code{p}=2
    corresponds to bivariate data.}
    \item{\code{k}}{The precise (simple) alternative hypothesis being tested.}
    \item{\code{cutoff1pct}}{A numeric vector consisting of the cut-off points for the \eqn{alpha}=0.01 level.}
    \item{\code{cutoff5pct}}{A numeric vector consisting of the cut-off points for the \eqn{alpha}=0.05 level.}
  }
}
\details{
  In order to test the significance of the Empirical Posterior Probability (EPP), it is
  necessary to generate its distribution under \eqn{H_0}, and compare the EPP to the sample quantile for the desired level
  of significance. However this simulation could take a considerably long time, and hence this
  table is provided to enable the experimenter to get a crude estimate of critical values at the
  \eqn{alpha}=0.05 and \eqn{alpha}=0.01 levels. The parameters under which the table was generated
  are part of the dataframe, allowing the experimenter to choose the closest set of conditions
  to his/her particular set-up.
}
\references{ 
  Fuentes, C. and Casella, G. (2009) Testing for the Existence of Clusters. 
  \emph{SORT} \bold{33}(2)
}
\examples{
data(cutoffs)
## maybe str(cutoffs) ; plot(cutoffs) ...
}
\keyword{datasets}
\keyword{ cluster }
\keyword{ htest }
\keyword{ multivariate }
