##########################################################################
# Version 5.0.0 of Clustering Algorithm
# Last edited Oct 09th 2008
##########################################################################
cluster.optimal <- function(data, nsim=1000, aR=0.4, p=2, k=2, a=2.01, b=0.990099, tau2=1, keep=4,
		mcs=0.2, file="", label="data") {

# Check arguments here ****************************************************
	if((k>=5)||(k<2)) stop("k can only take values 2, 3 or 4")
	if(mcs>=1) stop("mcs has to be a fraction between 0 and 1")
	if((a<=0)||(b<=0)||(tau2<=0)) stop("a, b and tau2 have to be non-negative")
	if(length(dim(data))!=2)
		stop("data has to be a 2D data array")
	if(dim(data)[2]!=p)
		stop("data columns do not conform to the given value of p")
# Arguments checked    ****************************************************

  n <- dim(data)[1]								# number of observations
  m <- round(mcs*n)								# minimum cluster size

  preConst <- p*a*log(2/b)-p*lgamma(a) 						# put this outside the k in 2:4 loop

  perm2 <- matrix(c(1,2,2,1), nrow=2, byrow=TRUE)
  perm3 <- matrix(c(1, 2, 3,
  1, 3, 2,
  2, 1, 3,
  2, 3, 1,
  3, 1, 2,
  3, 2, 1), nrow=6, byrow=TRUE)
  perm4 <- matrix(c(1, 2, 3, 4,
  1, 2, 4, 3,
  1, 3, 2, 4,
  1, 3, 4, 2,
  1, 4, 3, 2,
  1, 4, 2, 3,
  2, 1, 3, 4,
  2, 1, 4, 3,
  2, 3, 1, 4,
  2, 3, 4, 1,
  2, 4, 3, 1,
  2, 4, 1, 3,
  3, 2, 1, 4,
  3, 2, 4, 1,
  3, 1, 2, 4,
  3, 1, 4, 2,
  3, 4, 1, 2,
  3, 4, 2, 1,
  4, 2, 3, 1,
  4, 2, 1, 3,
  4, 3, 2, 1,
  4, 3, 1, 2,
  4, 1, 3, 2,
  4, 1, 2, 3), nrow=24, byrow=TRUE)
  perm <- list(perm2,perm3,perm4)

#**************************************************************************
  gdraw<-function(k) {
#**************************************************************************
    randind<-c(1,sample(c(rep(1,times=k-1),rep(0,times=(n-m*k)))))
    sind<-randind*c(1:(k+n-m*k))						#first indices of cluster start points
    sind<-sind[sind!=0]								#zeros eliminated
    start<-c(sind+c(0,1:(k-1))*(m-1),n)					#endpoint adjustment, end added
    temp<-rep(1,times=start[2]-start[1])
    for(i in 3:length(start))temp<-c(temp,rep(i-1,times=start[i]-start[i-1]))
    return(sample(c(temp,k)))		
  }
#**************************************************************************

#**************************************************************************
GroupMH <- function(cind, k) {
#**************************************************************************
means <- vector(mode="numeric", length=k*p)
variances <- vector(mode="numeric", length=k*p)
counts <- rep(0,4)
tmp <- .C("GroupMH", as.single(Y), as.integer(cind),as.integer(k), as.integer(p), 
	as.integer(length(cind)), as.single(means), as.single(variances), as.integer(counts), PACKAGE="bayesclust")
return(list(matrix(tmp[[6]], nrow=p, byrow=TRUE), matrix(tmp[[7]], nrow=p, byrow=TRUE), tmp[[8]][1:k]))
}
#**************************************************************************

#**************************************************************************
margM<-function(MV)
#**************************************************************************
{
N<-MV[[3]]
if (p>1) { for (i in 2:p) N <- rbind(N,MV[[3]]) }
return(sum((p*lgamma((MV[[3]]/2)+a))-(p/2)*log(MV[[3]]*tau2+1)) - sum(((N/2)+a)*log(N*MV[[2]]+(2/b))))
}
#**************************************************************************

#**************************************************************************
getRW<-function(cind) {
#**************************************************************************
mtest<-m
while(mtest<= m)
  {
  ind<-sample(n,1)
  mtest<-sum(cind==cind[ind])
  }
  temp<-c(1:k)[-cind[ind]]
  if(k>2) cind[ind]<-sample(temp,1)
  else cind[ind]<-temp
  return(cind)
}
#**************************************************************************

#**************************************************************************
getNew<-function(cind, logmarg1.old, Const) {
#**************************************************************************
# aR is the percentage of time it goes into a RW
# (1-aR) is the percentage of time it picks an independent draw
#**************************************************************************
if(runif(1)>aR)
  {
   Clustnew<-gdraw(k)
# use logmarg1.new <- margM(GroupMH(Clustnew, k))
   MVnew<-GroupMH(Clustnew, k)
  }  
else 
  {
   Clustnew <- getRW(cind)
   MVnew<-GroupMH(Clustnew, k)
  }
  lgnew <-lfact[k]-lchoose((n-(m*k)+k-1),k-1)-lfact.n +sum(lfactorial(MVnew[[3]]))
  lgold <-lfact[k]-lchoose((n-(m*k)+k-1),k-1)-lfact.n +sum(lfactorial(table(cind)))
#  MHR<-lgnew-lgold+log(Const + (1-aR)*exp(lgold))-log(Const + (1-aR)*exp(lgnew))
  logmarg1.new <- margM(MVnew)
  MHR <- logmarg1.new * (Const + (1 - aR)*exp(lgnew))^(-1) * (Const + (1 - aR)*exp(lgold)) * logmarg1.old^(-1)
  MHR <- min(1,MHR)
  if(runif(1)<MHR) return(list(Clustnew,logmarg1.new))
    else return(list(cind,logmarg1.old))
}
#**************************************************************************

#**************************************************************************
cluster_fix <- function(cluster, sizes, m) {
#**************************************************************************
takefrom <- which(sizes>m)
put.in <- which(sizes<m)
rsample <- NULL
for (i in takefrom) {
rsample <- c(sample(which(cluster==i), sizes[i]-m),rsample)
}
rsample <- sample(rsample)
for ( i in put.in ) {
put.in.i <- rsample[1:(m-sizes[i])]
cluster[put.in.i] <- i
rsample <- rsample[-(1:length(put.in.i))]
}
return(cluster)
}
#**************************************************************************

#**************************************************************************
cluster_relabel <- function(cluster.table, k) {
#**************************************************************************
top.clust <- cluster.table[1,1:n]
index <- array(TRUE, dim=c(n,k))

for (m in 2:keep) {
  tmp <- cluster.table[m,1:n]
  for (r in 1:k) {
    index[,r] <- tmp==r
  }

  max.cor <- -1 
  for (r in 1:gamma(k+1)) {
    for (t in 1:k) {
      tmp[index[,t]] <- perm[[k-1]][r,t] 
    }
    if (cor(top.clust,tmp) > max.cor) {
      max.cor <- cor(top.clust, tmp)
      cluster.table[m,1:n] <- tmp
    }
  }

}
return(cluster.table)
}
#**************************************************************************

  parameters <- list(n=n, k=k, min.clust.size=m, a=a, b=b, tau2=tau2, p=p, 
			random.walk=aR, iterations=nsim, dataset=data)
  
  Y <- scale(data)									#standardize

    logConst <- (k-1) * preConst          
#    set.seed(83.337 + k)

# The next few lines are purely to obtain the best cluster through kmeans.
# 'best' means smallest within SS
# That will serve as the starting point for the clustering algo.
    init.clust <- kmeans(Y, k)
    for (i in 2:10) {
      tmp <- kmeans(Y, k)    
      if (sum(tmp$withinss) < sum(init.clust$withinss))
        init.clust <- tmp
    }
    if (min(init.clust$size)>=m)
      Clust <- init.clust$cluster
    else 
      Clust <- cluster_fix(init.clust$cluster, init.clust$size, m)
    
# combine next two steps??
    MV1 <- GroupMH(Clust, k)
    logmarg1 <- margM(MV1)

    Const<-aR/(n*(k-1))

  if (keep>1) {
    Clusters <- array(0,dim=c(keep,n+1)) 
    Clusters[1,] <- c(Clust, logmarg1)

  lfact <- lfactorial(1:10)
  lfact.n <- lfactorial(n)

# populate the table with the first 'keep' clusters, making sure there are no duplicates 
    for(j in 2:keep) {
      temp <- getNew(Clust,logmarg1,Const)
      dup <- FALSE
      for(i in 1:(j-1)) {
        if (sum(Clusters[i,1:n]==temp[[1]])==n) {
          dup <- TRUE
          break
        }
      }
      # 'dup==FALSE' means that it is a new clustering/partition
      # 'dup==TRUE' means there is a duplicate clustering/partition already in the table
      if(dup==FALSE) 
        Clusters[j,] <- c(temp[[1]], temp[[2]])
      Clust <- temp[[1]]
      logmarg1 <- temp[[2]]
    }

# Store the minimum m(Y|w_k) and the lowest/smallest index at which it occurs
# Hence (*) will take the lowest/smallest index with the minimium logmarg1 value
    min.logmarg1 <- min(Clusters[,n+1])
    ind.min.logmarg1 <- min(which(Clusters[,n+1]==min.logmarg1))              # (*)

# to generate clusters according to MH and retain the best 'keep' clusters
    for(j in (keep+1):nsim) {
      temp<-getNew(Clust, logmarg1, Const)
      if (temp[[2]]>min.logmarg1) { 
        check.ind <- which(Clusters[,n+1]==temp[[2]])
        if(length(check.ind)==0) {
          Clusters[ind.min.logmarg1,] <- c(temp[[1]], temp[[2]])
          min.logmarg1 <- min(Clusters[,n+1])
          ind.min.logmarg1 <- min(which(Clusters[,n+1]==min.logmarg1))
        }
        else {
          # num.dup is the number of duplicate logmarg1's in the table
          num.dup <- length(check.ind)
          dup <- FALSE
          for(i in 1:num.dup) {
            if (sum(Clusters[check.ind[i],1:n]==temp[[1]])==n) {
              dup <- TRUE
              break
            }
          }
          # 'dup==FALSE' means that it is a new clustering/partition
          # 'dup==TRUE' means there is a duplicate clustering/partition already in the table
          if(dup==FALSE) {
            Clusters[ind.min.logmarg1,] <- c(temp[[1]], temp[[2]])
            min.logmarg1 <- min(Clusters[,n+1])
            ind.min.logmarg1 <- min(which(Clusters[,n+1]==min.logmarg1))
          }
        }
      }
      Clust <- temp[[1]]
      logmarg1 <- temp[[2]]
    }

   sort.idx <- sort(Clusters[,n+1], index=TRUE, decreasing=TRUE)$ix
   tmp2 <- cluster_relabel(Clusters[sort.idx,],k)
   Allclust <- list(clusters=tmp2[,1:n], logmarg=tmp2[,n+1])
  }
# end of if(keep>1) loop

# start of keep==1 loop
  else { 
    BestCluster <- c(Clust, logmarg1)
    for(j in 2:nsim) {
      temp<-getNew(Clust, logmarg1, Const)
      if (temp[[2]] > BestCluster[n+1]) {
        BestCluster <- c(temp[[1]], temp[[2]])
      }
      Clust <- temp[[1]]
      logmarg1 <- temp[[2]]
    }
   tmp2 <- BestCluster
   Allclust <- list(cluster=tmp2[1:n], logmarg=tmp2[n+1])
  }

      lbl <- paste(label,1,sep="")
      cluster.optimal.obj <- list(param=parameters, Allclust) 
      names(cluster.optimal.obj)[2] <- lbl

  class(cluster.optimal.obj) <- "cluster.optimal"

  if (nchar(file)!=0) {
    save(cluster.optimal.obj, file=file) 
  }
  cluster.optimal.obj

}

