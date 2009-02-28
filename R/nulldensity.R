##########################################################################
# Version 3.0.0
# Last edited Aug 30th 2008
##########################################################################
# nsim	: number of simulations, i.e., number of posterior probabilities to generate
# n 	: number of observations in dataset
# k	: the hypothesis to be carried out - either 2,3 or 4
# mcs	: minimum cluster size ( a fraction from 0 to 1)
# a,b, tau2	: prior parameters, all positive
# sig2	: prior parameter, should be a vector of length p, all positive
# prop	: proportion of partition-space to sample ( a fraction from 0 to 1)
# p	: each observation will be a p-variate vector (number of PCAs)
# file	: filename to save created R object to. The R object is saved,
# 	  not a txt file.
##########################################################################
nulldensity <- function(nsim, n, k, mcs=0.2, a=2.01, b=0.990099, 
	tau2=1, prop=0.25, p, file=""){

# Check arguments here ****************************************************
	if((mcs>=1)||(mcs >= 1/k)) stop("mcs has to be a fraction between 0 and 1/k")
	if((a<=0)||(b<=0)||(tau2<=0)) stop("a, b and tau2 have to be non-negative")
	sig2=rep(1,p) 
#	if((length(sig2)!=p)||(min(sig2)<=0)) 
#		stop("sig2 has to be a vector of positive values, of length p")
# Arguments checked    ****************************************************

#**************************************************************************
num.part <- function(n,k) {
#**************************************************************************
	M <- matrix(0,nrow=n+1,ncol=n+1)
	M[1,1] <- 1
	for(i in 1:n){ 	
	 for(j in 1:i){M[i+1,j+1] <- M[i,j]+M[i+1-j,j+1]}
         }
          M[n+1,k+1]
       }
#**************************************************************************

#**************************************************************************
loggammak1 <- function(n, k, a, b, sig2, tau2, p) {
#**************************************************************************
	logconst <- ((p*a*(k-1))*log(2/b) + (p/2)*log(n*tau2+1) - p*(k-1)*lgamma(a) - p*lgamma(n/2 + a))
	return (logconst)
	}
#**************************************************************************

#**************************************************************************
loggammak2 <- function(n, N, k, a, tau2, p) {
#**************************************************************************
	logA <- p*sum(lgamma(N/2 + a)) - (p/2)*(sum(log(N*tau2 + 1)))
	return(logA)
	}
#**************************************************************************

#**************************************************************************
loggammak3 <- function(n, obs, N, k, a, b, sig2, tau2) {
#**************************************************************************
	logB <- ((1-k)*a)*log(sig2)
	logC <- (n/2 + a)*log(sum(obs[-n]) + 2/(b*sig2))
	nul <- which(N==1)
	t <- length(nul)
        sums <- rep(0, k)
	if(t==0) {
		index <- rep(1:k, times=N-1)
		logD <- sum((N/2+a)*
			log(.C("nullsum", as.single(obs[1:(n-k)]), as.integer(index), as.integer(length(index)), 
			as.integer(k), as.single(sums), PACKAGE="bayesclust")[[5]][1:k] + 2/(b*sig2)))
		}
	else {
		index <- c(rep(1:k, times=N-1),nul)
		 tmp1 <- c(obs[1:(n-k)],rep(0,t))	
		 logD <- sum((N/2+a)*
			 log(.C("nullsum", as.single(tmp1), as.integer(index), as.integer(length(index)), 
			 as.integer(k), as.single(sums), PACKAGE="bayesclust")[[5]][1:k] + 2/(b*sig2)))
		}

	return (logB + logC - logD)
	}
#**************************************************************************

#**************************************************************************
bayesfactor.sampledM <- function(n, obs, k, m, a, b, sig2, tau2, sample.num, p) {
#**************************************************************************

	# All the small gamma values will be stored in this BF array
	BF <- vector(mode="numeric", length=sample.num)


	#begin loop

	for (j in 1:sample.num) {
		randind <- c(1,sample(c(rep(1,times=k-1),rep(0,times=(n-m*k)))))
		sind <- randind*c(1:(k+n-m*k))              		#first indices of cluster start points
		sind <- sind[sind!=0]                           	#zeros eliminated
		start <- c(sind+c(0,1:(k-1))*(m-1),n+1)          	#endpoint adjustment, end added
		lengthstart <- length(start)
		N <- rep(0,times=lengthstart-1)
		for(u in 1:(lengthstart-1)) N[u] <- start[u+1]-start[u]
                logterm2 <- loggammak2(n, N, k, a, tau2, p)
	        logterm3 <- vector("numeric",p)
		if (p>1) 
                  { for (r in 1:p) logterm3[r] <- loggammak3(n, obs[r,], N, k, a, b, sig2[r], tau2)}	# key line!!!
                else
                  { logterm3 <- loggammak3(n, obs, N, k, a, b, sig2, tau2)} 
	        BF[j] <- exp(logterm1 + logterm2 + sum(logterm3))
	}

	p_H0 <- 1/(1+ mean(BF))
	return(p_H0)
	}
#**************************************************************************

	Y <- matrix(0, p*nsim, n)
	  for(i in 1:(p*nsim)) Y[i,] <- rchisq(n,1)
	
	sample.num <- round(prop * num.part(n,k))
	m <- max(round(mcs * n),1)

	Pvalue <- vector(mode="numeric", length=nsim)

	# Compute the constant T
	logterm1 <- loggammak1(n, k, a, b, sig2, tau2, p)

	for(i in 1:nsim) Pvalue[i] <- bayesfactor.sampledM(n, Y[((i-1)*p+1):(i*p),], k, m, a, b, sig2, tau2, sample.num, p)
	orderstat <- sort(Pvalue)

	parameters <- list(n=n, k=k, min.clust.size=m, a=a, b=b, sig2=sig2, tau2=tau2, prop.sampled=prop, p=p)
	nulld <- list(param=parameters, gen.values=orderstat)
	
	class(nulld) <- "nulldensity"

	if (nchar(file)!=0) {
	  save(nulld, file=file) 
	}

	nulld

}
