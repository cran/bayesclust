print.emp2pval <- function(x, ...) {
  if(class(x)!="emp2pval")
    stop("The first argument is not of the right class. Please see the help pages for more information")

  print(x$pvals)
}

