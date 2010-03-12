summary.cluster.test.reps <- function(object, ...) {
  x <- object[[1]]
  data.used <- names(x)[3]
  iterations <- tail(x$iterations, n=1) 
  num.obs <- x$param$n
  m <- x$param$min.clust.size
  p <- x$param$p
  null <- x$param$k

  rownames <- paste("rep", 1:length(object), sep="")
  postprobs <- vector("numeric", length=length(object))
  for (r in 1:length(object) ) {
    postprobs[r] <- tail(object[[r]][[3]], n=1)
  }
  TAB <- data.frame("Post.Probs"=round(postprobs, 4), row.names=rownames)
  res <- list( data.used=data.used , iterations=iterations,  num.obs=num.obs,  
    m=m, p=p, null=null, probs=TAB, replications=length(object))
  class(res) <- "summary.cluster.test.reps"
  res
}
