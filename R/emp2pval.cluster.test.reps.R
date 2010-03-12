emp2pval.cluster.test.reps <- function(x, y, ignore=FALSE) {
  if(class(y)!="nulldensity")
    stop("The second argument is not of the right class. Please see the help 
      pages for more information")

  tmp01 <- x[[1]]
  if(ignore==FALSE) {
    if((tmp01$param$n!=y$param$n)||(tmp01$param$k!=y$param$k)||
      (tmp01$param$min.clust.size!=y$param$min.clust.size)||
      (tmp01$param$a!=y$param$a)||(tmp01$param$b!=y$param$b)||
      (tmp01$param$tau2!=y$param$tau2)||(tmp01$param$p!=y$param$p))
      stop("parameters of cluster test do not match parameters used when 
        generating null density!")
  }


  row.names <- names(tmp01)[3]

  mean.emp.prob <- 0
  for (r in 1:length(x)) {
    mean.emp.prob <- mean.emp.prob + tail(x[[r]][[3]], n=1)
  } 
  mean.emp.prob <- mean.emp.prob/length(x)

  pvalue <- findInterval(mean.emp.prob, y$gen.values) + 1
  pvalue <- pvalue/length(y$gen.values)  
  pvalue <- min(1, pvalue)
  
  pvals <- data.frame(mean.emp.prob, pvalue, row.names=row.names)
  res <- list(param=tmp01$param, pvals=pvals)
  class(res) <- "emp2pval"
  res
}
