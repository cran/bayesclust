plot.cluster.optimal <- function(x, varToPlot=c(1,2), clustToDisp=c(1,2,3,4), 
  ...) {
# checking number of dimensions and which of them are to be plotted.
  p <- x$param$p
  if (length(varToPlot)!=2) { 
    stop("argument varToPlot should be a vector of length 2, indicating which 
      of the 'p' variables to plot")
  }
  if (max(varToPlot)>p) {
    stop(paste("varToPlot argument is asking to plot variable", max(varToPlot), 
      "whereas there are only", p, "variables in the data.", sep=" "))
  }

# checking how many clusterings were kept and which are to be plotted.
  data <- x$param$dataset
  if (length(dim(x[[2]][[1]]))==0) keep <- 1
  else keep <- dim(x[[2]][[1]])[1] 
  if (length(clustToDisp)>4) {
    warning("This function only displays 4 clusterings at a time. Please call 
      the function again, in order to plot the remaining clusters.")
    clustToDisp <- clustToDisp[1:4]
  }
  if (max(clustToDisp)>keep) {
    stop(paste("Sorry, only", keep, "clusterings were kept.", sep=" "))
  }

  if(keep==1) 
    plot(data[,varToPlot], col=x[[2]][[1]], main="Optimal Cluster", 
      xlab=paste("x",varToPlot[1],sep=""), 
      ylab=paste("x",varToPlot[2],sep=""), ...)
  else {
    par(mfrow=c(2,2))
    for (i in 1:min(4, length(clustToDisp))) {
      mainTitle <- paste("Optimal Clusters\nRank", clustToDisp[i], sep=" ")
      plot(data[,varToPlot], col=x[[2]][[1]][clustToDisp[i],], 
        main=mainTitle, xlab=paste("x",varToPlot[1],sep=""), 
        ylab=paste("x",varToPlot[2],sep=""), ...)
    }
  }
  
}
