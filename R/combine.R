combine <- function(... ) {
  argList <- match.call()
  argList <- as.list(argList)
  argList[[1]] <- NULL
  numModels <- length(argList)
  if (numModels==1) 
    stop("At least 2 tests have to be run in order to combine")

  objList <- vector("list", length=numModels)
  for (i in 1:numModels ) {
    objList[[i]] <- get(as.character(argList[[i]]))
    if((class(objList[[i]])!="cluster.test")&&
      (class(objList[[i]])!='cluster.test.reps'))
        stop(paste("input arguments have to be objects returned from ", 
          "function 'cluster.test()'", sep=""))
  }

# check if all parameters are same for all objects
  tmpList <- vector("list", length=numModels)
  for (i in 1:numModels ) {
    if (class(objList[[i]])=="cluster.test") {
      tmpList[[i]] <- objList[[i]]$param
      tmpList[[i]][[2]] <- NULL
    }
    else  {
      tmpList[[i]] <- objList[[i]][[1]]$param
      tmpList[[i]][[2]] <- NULL
    }
  }
  for (i in 2:numModels ) {
    if (!identical(tmpList[[1]], tmpList[[i]])) {
      stop("the objects should have been generated with the same parameters")
    }
  }

  hypTestedK <- vector("numeric", length=numModels)
  bayesFactorK1 <- vector("numeric", length=numModels)
  for (i in 1:numModels ) {
    if (class(objList[[i]])=="cluster.test") {
      hypTestedK[i] <- objList[[i]]$param[[2]]
      bayesFactorK1[i] <- (1/tail(objList[[i]]$data, n=1)) - 1
    }
    else {
      hypTestedK[i] <- objList[[i]][[1]]$param[[2]]
      tmpSum <- 0
      numReps <- length(objList[[i]])
      for (j in 1:numReps ) {
        tmpSum <- tmpSum + (1/tail(objList[[i]][[j]]$data, n=1)) - 1
      }
      bayesFactorK1[i] <- tmpSum/numReps
    }
  }
  if (length(unique(hypTestedK))<length(hypTestedK))
    stop("Duplicated hypothesis tests have been given")

  hypTestedK <- c(1, hypTestedK)
  bayesFactorK1 <- c(1, bayesFactorK1)
  postProb <- bayesFactorK1/sum(bayesFactorK1)

  rownames <- paste("K=", hypTestedK, sep="")
  data.frame("Postr.Prob"=round(postProb, digits=8), row.names=rownames)
}
