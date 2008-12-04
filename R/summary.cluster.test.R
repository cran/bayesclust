summary.cluster.test <- function(object, ...) {
  data.used <- names(object)[3]
  iterations <- tail(object$iterations, n=1) 
  num.obs <- object$param$n
  m <- object$param$min.clust.size
  p <- object$param$p
  null <- object$param$k

    rownames <- names(object[3])
    postprobs <- tail(object[[3]], n=1)
  TAB <- data.frame("Post.Probs"=round(postprobs, 4), row.names=rownames)
  res <- list( data.used=data.used , iterations=iterations,  num.obs=num.obs,  m=m, p=p, null=null, probs=TAB)
  class(res) <- "summary.cluster.test"
  res
}
