
DP_fn.smart_split <- function(parm, data)
{	
	# pick a large cluster
	# singletons never picked
	pick.g <- sample(1:parm$clust$G, size=1, prob=(parm$clust$C.m.vec-1)^(2))

	# which members of this group look very different from the "center"
	# i.e. have small correlations (-ve is worse than close to 0)
	indx.g <- which(parm$clust$c.v==pick.g)
	cor.v <- as.vector(cor(parm$X[,indx.g], parm$clust$A.mt[,pick.g]))
	
	##diff.indx.g <- indx.g[which.min(cor.v)]

	# pick "num.picked" members with smallest correlations less than min.cor
	num.picked <- 2
	min.cor <- .95
	cutoff <- quantile(cor.v, prob=min(1,5/length(cor.v)))
	cutoff <- min(min.cor, cutoff)
	diff.indx.g <- indx.g[which(cor.v < cutoff)]

	# try to split group members with R wrt center less than min.cor
	##min.cor <- .5
	##diff.indx.g <- indx.g[which(cor.v < min.cor)]

	parm$split.changed.v <- NULL
	
	if (length(diff.indx.g)>0)
		{
	 	for (cc in diff.indx.g)
			{previous.parm <- parm

			parm$k <- cc

			tmp <- DP_fn.gibbs(k=parm$k, parm, data, no.prior=FALSE)
			parm <- tmp[[1]]	
			}	

		parm$split.changed.v <- parm$clust$c.v[diff.indx.g] != pick.g
		}

 	err <- DP_fn.consistency.check(parm)
	if (err > 0)
		{stop(paste("END SPLIT: failed consistency check: err=",err))
		}

	parm
}


###########################################################

DP_fn.log.lik <- function(...) {
    stop("Missing function")
}

DP_fn.consistency.check <- function(...) {
    stop("Missing function")
}

DP_fn.gibbs <- function(...) {
    stop("Missing function")
}


DP_fn.smart_merge <- function(parm, data)
{	

	old.parm <- parm

	# pick two small clusters
	pick.g <- sample(1:parm$clust$G, size=2, prob=parm$clust$C.m.vec^(-1))
	indx.g <- which(parm$clust$c.v %in% pick.g)
	old.c.v <- parm$clust$c.v[indx.g]
	ratio.cor2.v <- as.vector(cor(parm$X[,indx.g], parm$clust$A.mt[,pick.g[1]])^2)/as.vector(cor(parm$X[,indx.g], parm$clust$A.mt[,pick.g[2]])^2)
	prob.v <- ratio.cor2.v/(1+ratio.cor2.v)

	new.c.v <- NULL
	for (xxx in 1:length(indx.g))
		{new.c.v <- c(new.c.v, sample(pick.g, size=1, prob=c(ratio.cor2.v[xxx],1)))
		}

	switched.flag <- old.c.v!=new.c.v
	log.prop <- 0

	parm$clust$c.v[indx.g] <- new.c.v

	for (gg in pick.g)
			{flag.gg <- (new.c.v==gg)
			parm$clust$C.m.vec[gg] <- sum(flag.gg)
			}

	parm$merge.flip <- -1

	if (sum(switched.flag) > 0)
		{log.prop <- log.prop + sum(log(prob.v[switched.flag])) 

		I.k <- indx.g[switched.flag]
		short.X <- matrix(parm$X[,I.k], ncol=length(I.k))

		new.indx <- which(parm$clust$C.m.vec>0)
		old.indx <- which(old.parm$clust$C.m.vec>0)

		new.G <- length(new.indx)
		old.G <- length(old.indx)

		# copied from clust,functions.R based on Lijoi and Prunster
		new.log.lik <- sum(log(parm$b1 + (1:(new.G-1))*parm$d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-parm$d), (parm$clust$C.m.vec[new.indx]-1)))
		old.log.lik <- sum(log(parm$b1 + (1:(old.G-1))*parm$d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-parm$d), (old.parm$clust$C.m.vec[old.indx]-1)))

		log.true <- as.numeric((parm$clust$C.m0>0)-(old.parm$clust$C.m0>0))*log(parm$b0) + new.log.lik - old.log.lik
		
		### 
		new.log.lik <- 0
		short.c.v <- new.c.v[switched.flag]
		#
		for (gg in sort(unique(short.c.v)))
			{flag.gg <- (short.c.v==gg)
			if (sum(flag.gg)>0)
				{x.mt <- matrix(short.X[,flag.gg], ncol=sum(flag.gg))
				new.log.lik <- new.log.lik + sum(DP_fn.log.lik(gg, x.mt, parm, colSums=TRUE))
				}
			}

		### 

		old.log.lik <- 0
		short.c.v <- old.c.v[switched.flag]
		#
		for (gg in sort(unique(short.c.v)))
			{flag.gg <- (short.c.v==gg)
			if (sum(flag.gg)>0)
				{x.mt <- matrix(short.X[,flag.gg], ncol=sum(flag.gg))
				old.log.lik <- old.log.lik + sum(DP_fn.log.lik(gg, x.mt, old.parm, colSums=TRUE))
				}
			}

		###

		log.true <- log.true + new.log.lik - old.log.lik

			
		################
		## T O S S
		################

		prob <- exp(min((log.true-log.prop),0))
		flip <- as.logical(rbinom(n=1, size=1, prob=prob))
		if (!flip) {parm <- old.parm}

 		err <- DP_fn.consistency.check(parm)
		if (err > 0)
			{stop(paste("END STEP: failed consistency check: err=",err))
			}

		parm$merge.flip <- flip

		} # end massive if (sum(switched.flag) > 0)

	parm
}



###########################################################


DP_fn.smart_split_merge <- function(parm, data)
{	
	parm$merge.flip <- NA

	if (parm$clust$G < parm$G.max)
		{parm <- DP_fn.smart_split(parm, data)
		}

	parm <- DP_fn.smart_merge(parm, data)

	parm
}


