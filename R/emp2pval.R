emp2pval <- function(x, y, ignore=FALSE) {
  if((class(x)!="cluster.test")||(class(y)!="nulldensity"))
    stop("The first two arguments are not of the right class. Please see the help pages for more information")

  if(ignore==FALSE) {
    if((x$param$n!=y$param$n)||(x$param$k!=y$param$k)||(x$param$min.clust.size!=y$param$min.clust.size)||
	(x$param$a!=y$param$a)||(x$param$b!=y$param$b)||(x$param$tau2!=y$param$tau2)||(x$param$p!=y$param$p))
      stop("parameters of cluster test do not match parameters used when generating null density!")
  }


  row.names <- names(x)[3]

    emp.prob <- tail(x[[3]], n=1)

  pvalue <- findInterval(emp.prob, y$gen.values) + 1
  pvalue <- pvalue/length(y$gen.values)  
  pvalue <- min(1, pvalue)
  
  pvals <- data.frame(emp.prob, pvalue, row.names=row.names)
  res <- list(param=x$param, pvals=pvals)
  class(res) <- "emp2pval"
  res
}
