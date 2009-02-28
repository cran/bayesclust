fdr.test <- function(namelist, q=0.05) {
  
  len.namelist <- length(namelist)

  for (x in namelist) {
    tmp1 <- get(x)
    if(class(tmp1)!="emp2pval")
      stop("Arguments should be of class emp2pval. Please see the help pages for more information")
  }
  if ((q<=0)||(q>1))
      stop("q should be between 0 and 1. Please see the help pages for more information")

  fdr.seq <- q * ((1:len.namelist)/len.namelist)
  
  k <- NULL
  labels <- NULL
  emp.pvals <- NULL
  for (x in namelist) {
    tmp1 <- get(x)
    k <- c(k, tmp1[[1]][[2]])
    labels <- c(labels, x)
    emp.pvals <- rbind(emp.pvals, tmp1[[2]])
  }

  consol <- data.frame(labels, k, emp.probs=round(emp.pvals[,1],4), pvalues=round(emp.pvals[,2],4))
  
  ord.ind <- order(consol$pvalues)
  tmp <- which(sort(consol$pvalues)<=fdr.seq)
  if(length(tmp)==0) {
      cat("None of the hypothesis are significant when controlling FDR at level alpha = ", q ,"\n")
  }
  else {
  sig.hyp <- max(which(sort(consol$pvalues)<=fdr.seq))
    if(sig.hyp==len.namelist) {
      cat("All the hypotheses are significant when controlling FDR at level alpha = ", q ,"\n")
    }
    else {
      cat("Significant clusters were found in the following datasets when controlling FDR at level alpha = ", q ,"\n")
      sig.clust <- consol[ord.ind[1:sig.hyp], ]
      print(sig.clust)
    }
  }
}
