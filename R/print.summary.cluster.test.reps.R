print.summary.cluster.test.reps <- function(x, ...) {
  cat("\n")
  cat("Cluster test conducted on data object ", x$data.used, ", with ", x$iterations, 
    " iterations.\n", sep="")
  cat("Num. observations\t: ", x$num.obs, "\n", sep="")
  cat("Min cluster size\t: ",x$m,"\n", sep="")
  cat("p\t\t\t: ", x$p,"\n",sep="")
  cat("H1\t\t\t: k = ",x$null,"\n",sep="")

  cat("**************************************** \n")
  cat("Final Empirical Posterior Probabilities: \n")
  cat("**************************************** \n")
  print(x$probs)
  cat("\n")
  cat("Please run emp2pval to obtain the corresponding P-values for the mean 
    of the above statistics\n")
}
