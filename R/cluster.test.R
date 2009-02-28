##########################################################################
# Version 8.0.0
# Last edited Oct 09th 2008
##########################################################################
cluster.test <- function(data, nsim=1000, aR=0.4, p=2, k=2, a=2.01, b=0.990099, tau2=1, 
		mcs=0.2, file="", label="data") {

# Check arguments here ****************************************************
	if((mcs>=1)||(mcs >= 1/k)) stop("mcs has to be a fraction between 0 and 1/k")
	if((a<=0)||(b<=0)||(tau2<=0)) stop("a, b and tau2 have to be non-negative")
	if(length(dim(data))!=2)
		stop("data has to be in the form of a 2D data array")
	if(dim(data)[2]!=p)
		stop("data columns do not conform to the given value of p")
# Arguments checked    ****************************************************

  n <- dim(data)[1]								# number of observations
  m <- max(round(mcs*n),1)							# minimum cluster size

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
counts <- rep(0,k)
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
# This part essentially uses the Fundamental Theorem of Simulation,
# because it samples uniformly from the bigger set, but only keeps 
# those from the smaller one (i.e., those clusters that can spare the
# extra observation.
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
#**************************************************************************
getNew<-function(cind,MV,Const,lgold) {
#**************************************************************************
# aR is proportion of time the algo goes into RW. Make it high!!
if(runif(1)>aR)
  {
   Clustnew<-gdraw(k)
   MVnew<-GroupMH(Clustnew, k)
  }  
else 
  {
   Clustnew <- getRW(cind)
   MVnew<-GroupMH(Clustnew, k)
  }
  counts.old <- as.vector(table(cind))
  n.old <- sum(counts.old[counts.old>m])
  counts.new <- as.vector(table(Clustnew))
  n.new <- sum(counts.new[counts.new>m])
  lgnew<-lfact.k-lchoose((n-(m*k)+k-1),k-1)-lfact.n +sum(lfactorial(MVnew[[3]]))
  MHR<-lgnew-lgold+log(Const/n.new + (1-aR)*exp(lgold))-log(Const/n.old + (1-aR)*exp(lgnew))
  MHR<-min(1,exp(MHR))
  if(runif(1)<MHR) return(list(Clustnew,lgnew,MVnew))
    else return(list(cind,lgold,MV))
}
#**************************************************************************

#**************************************************************************

  parameters <- list(n=n, k=k, min.clust.size=m, a=a, b=b, tau2=tau2, p=p, data=data)
  iter <- seq(from=nsim%%500, to=nsim, by=500)
  if (iter[1]==0) iter <- iter[-1]
  
  Y <- scale(data)								#standardize
  ClusterStat <- vector("numeric", nsim) 
  preConst <- p*a*log(2/b)-p*lgamma(a) 					# put this outside the k in 2:4 loop
  logmarg0 <- p*lgamma(n/2 + a)-(p/2)*log(n*tau2+1) - (n/2 + a)*sum(log(n*diag(var(Y))+ 2/b))
  T1 <- vector("numeric", nsim)
  lfact.k <- lfactorial(k)
  lfact.n <- lfactorial(n)

    logConst <- (k-1) * preConst          

# combine next two steps??
    Clust <- gdraw(k)
    MV1 <- GroupMH(Clust, k)
    logmarg1 <- margM(MV1)
    T1[1] <- logmarg1

#    Const<-aR/(n*(k-1))
    Const<-aR/(k-1)
    log.g<-lfact.k-lchoose((n-(m*k)+k-1),k-1)-lfact.n +sum(lfactorial(Clust))
    
    for(j in 2:nsim) {
      temp <- getNew(Clust, MV1, Const, log.g)
      MV1 <- temp[[3]]
      T1[j] <- margM(MV1) #(=logmarg1)
      log.g<-temp[[2]]
      Clust <- temp[[1]]
    }

    topG <- cumsum( exp(T1 + logConst - logmarg0))
    BFG <- topG/c(1:nsim)
    ClusterStat <- 1/(1+BFG) 
    ClusterStat <- ClusterStat[iter]

      lbl <- paste(label,1,sep="")
      cluster.test.obj <- list(param=parameters, iterations=iter, ClusterStat) 
      names(cluster.test.obj)[3] <- lbl

  class(cluster.test.obj) <- "cluster.test"

  if (nchar(file)!=0) {
    save(cluster.test.obj, file=file) 
  }
  cluster.test.obj
}

