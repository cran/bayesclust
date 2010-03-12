plot.cluster.test.reps <- function(x, ...) {
  numReps <- length(x)

  if (numReps <= 4) {
    par(mfrow=c(2,2))
    for (r in 1:numReps ) {
      y <- x[[r]]
      plot(y$iterations, y[[3]], ylim=c(0,1), type="l", xlab=names(y)[2], 
        ylab="Running Posterior Prob", 
        main=paste("Replication number",r , sep=" "), ...)
    }
  }
  else {
    par(mfrow=c(2,2))
    for (r in 1:4 ) {
      y <- x[[r]]
      plot(y$iterations, y[[3]], ylim=c(0,1), type="l", xlab=names(y)[2], 
        ylab="Running Posterior Prob",
        main=paste("Replication number",r , sep=" "), ...)
    }
    devAskNewPage(ask=TRUE)
    for (r in 5:numReps ) {
      y <- x[[r]]
      plot(y$iterations, y[[3]], ylim=c(0,1), type="l", xlab=names(y)[2], 
        ylab="Running Posterior Prob",
        main=paste("Replication number",r , sep=" "), ...)
      if(r%%4==0) {
        devAskNewPage(ask=TRUE)
      }
    }
  }

}
