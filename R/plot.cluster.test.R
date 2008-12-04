plot.cluster.test <- function(x, ...) {
  if(class(x)!="cluster.test")
    stop("The first argument is not of the right class. Please see the help pages for more information")

    plot(x$iterations, x[[3]], ylim=c(0,1), type="l", xlab=names(x)[3], ylab="Running Posterior Prob", ...)
}
