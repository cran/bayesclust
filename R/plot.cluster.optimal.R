plot.cluster.optimal <- function(x, ...) {
  p <- x$param$p
  if ((p<2)||(p>6))
    stop("This plot function only works for p between 2 and 6.")

  data <- x$param$dataset
  if (length(dim(x[[2]][[1]]))==0) keep <- 1
  else keep <- dim(x[[2]][[1]])[1] 
  
  if(p==2) {
      if(keep==1) 
      plot(data[,1:2], pch=20, col=x[[2]][[1]], 
	main="Optimal Cluster", xlab="x1", ylab="x2", ...)
      else 
      plot(data[,1:2], pch=20, col=x[[2]][[1]][1,], 
	main="Optimal Cluster", xlab="x1", ylab="x2", ...)
  }
  else {
      if (p==3) {
        par(mfrow=c(3,1))
      }
      else if (p==4) {
	par(mfrow=c(2,3))
      }
      else if (p==5) {
	par(mfrow=c(2,5))
      }
      else {
	par(mfrow=c(3,5))
      }
      if(keep==1) {
       for (j in 1:p) 
	for (k in (j+1):p) {
	  xvar <- j
	  yvar <- k
          plot(data[,c(xvar, yvar)], pch=20, col=x[[2]][[1]], 
	    main="Optimal Cluster", 
	      xlab=paste("x",xvar, sep=""), ylab=paste("x",yvar,sep=""), ...)
       }
      }
      else {
       for (j in 1:(p-1)) 
	for (k in (j+1):p) {
	  xvar <- j
	  yvar <- k
          plot(data[,c(xvar, yvar)], pch=20, col=x[[2]][[1]][1,], 
	    main="Optimal Cluster", xlab=paste("x",j, sep=""), ylab=paste("x",k,sep=""), ...)
       }
      }
  }
  
}
