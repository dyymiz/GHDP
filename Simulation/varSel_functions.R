
################################
### explanation of the different X matrices
##################################

## Expand and knots calculation in spline.
## Matrix orthogonalization, Gram-Schmidt algorithm.
## Calculate likelihoods  L1.v and L3.v

## n2 is the same as n in clustering part. Here, n2=n1+n0, n0=number of training cases

## parm$all.X is the same as All.Stuff$parm$X. It is full n2 X p matrix of actual/imputed covariates

## parm$X: DOES NOT EXIST for responses analysis. It is replaced by parm$all.X

## parm$clust$X.mt: n2 by (parm$clust$G+1) matrix of one covariate from each of the parm$clust$G 
## clusters plus intercept. It is computed in fn.pick.predictors. It is indep of num.knots (num.knots?)

## parm$X.mt is the collection of INTERCEPT followed by parametric PREDICTORS stratified by knots 
## followed by SECOND ORDER of parametric PREDICTORS stratified by knots (order of types extremely 
## important): It is n2 by (1+(parm$num.knots+1)*parm$clust$rho+(parm$num.knots+3)*parm$clust$rho_3) 
## matrix. Computed in fn.get.X.mt

## Matching the columns of parm$X.mt[,-1], parm$clust.beta.v is now the regression coefficients of 
## parametric PREDICTORS stratified by knots followed by SECOND ORDER of parametric PREDICTORS 
## stratified by knots (order of types extremely important): 
## It is vector of length (parm$num.knots+1)*parm$clust$rho+(parm$num.knots+1)*parm$clust$rho_2

## parm$clust$gamma.v is indep of num.knots; so is parm$clust$rho which is sum(parm$clust$gamma.v)

## parm$clust$gamma_2.v is indep of num.knots; so is parm$clust$rho_2 which is sum(parm$clust$gamma_2.v)
## parm$clust$rho: number of linear predictors
## parm$clust$rho_2: number of second order predictors
## parm$clust$rho_3: number of second order predictors

######################################



fn_var.sel_eda <- function(parm, data)
{	
  parm = fn_gen.Y.summary(parm)
  
  tmp <-  aggregate(t(parm$all.X) ~ parm$clust$c.v, FUN=mean, simplify=TRUE)
  parm$clust$A.mt <- as.matrix(t(tmp[,-1]))
  dimnames(parm$clust$A.mt) = NULL
    
 
	########################
 
	parm <- fn.init.pick.predictors(parm)

	################

	parm
}


fn.mboost.init <- function(parm, data)
	{	
	# code copied from "mboost.R"
	# with some modifications

	options(warn=0)
	library(mboost)
	options(warn=2)

	ind.tr <- data$non.missing.indx
	ind.test <- setdiff(1:n2, data$non.missing.indx)

	df <- data.frame(exp(data$Y), data$delta, parm$all.X[ind.tr,parm$small.indx])

	names(df)[1:2] <- c("t", "di")

	t <- df$t
	di <- df$di
	x <- df[,-(1:2)]

	options(warn=0)
	train.obj <- glmboost(Surv(t, di) ~ ., data=df, family = CoxPH(), center=TRUE)
	options(warn=2)

	vec <- names(coef(train.obj))[-1]
	selected.cols <- as.numeric(unlist(strsplit(vec, split="X"))[seq(2,by=2,len=length(vec))])

	selected.cols
	}




fn_varsel.init <- function(true, data, All.Stuff, num.knots){

	parm <- All.Stuff$parm
	
	############################################
	# parm$X is renamed parm$all.X in this part of the procedure
	############################################
	parm$all.X <- parm$X
	parm$X <- NA

	parm <- fn_var.sel_eda(parm, data)

	##############

	parm$prior$gamma.g$omega <- parm$prior$gamma.g$omega_3 <-.1 

	############################################################
	# Option 1: using boosting and parm$small.indx (1 from each cluster)
	# 		to get predictors
	############################################################

	test.flag <- try(fn.mboost.init(parm, data))

	# survival regression using only selected cols
	if (mode(test.flag) == "numeric")
		{
		parm$selected.cols <- selected.cols <- test.flag
	 	tmp <- survreg(Surv(exp(data$Y), data$delta) ~ parm$all.X[data$non.missing.indx,parm$small.indx[selected.cols]], control = list(maxiter=20000), dist="gaussian")
	 	tmp2 <- summary(tmp) 
		table <- tmp2$table[-(length(selected.cols)+2),]
		### plot(tmp$linear.pred, data$Y); 

		parm$clust$gamma.v <- parm$clust$beta.v <- array(0, parm$clust$G)
		parm$clust$gamma.v[selected.cols] <- 1
		parm$clust$beta.v[selected.cols] <- as.numeric(table[-1,1])
		parm$clust$beta0 <- as.numeric(table[1,1])
		}
	


	############################################################
	# Option 2: Boosting (Option 1) fails and parm$clust$G >= n0
	############################################################
	# single-variable survival regression for each candidate column
	# Pick "init.rho" clusters having lowest separate p-values
	# survreg.control(maxiter=20000)
	
	if (mode(test.flag) != "numeric")
		{
		p.val.v <- array(,parm$clust$G)

		for (ii in 1:parm$clust$G)
			{tmp <- survreg(Surv(exp(data$Y), data$delta) ~ parm$all.X[data$non.missing.indx,parm$small.indx[ii]], control = list(maxiter=20000), dist="gaussian")
			tmp2 <- summary(tmp) 
			p.val.v[ii] <- tmp2$table[2,4]
			}
		init.rho <- round(min(c(length(data$non.missing.indx), .1*parm$clust$G)))
		cutoff <- quantile(p.val.v, prob=init.rho/parm$clust$G)
		selected.cols <- which(p.val.v <= cutoff)

		tmp <- survreg(Surv(exp(data$Y), data$delta) ~ parm$all.X[data$non.missing.indx,parm$small.indx[selected.cols]], control = list(maxiter=20000), dist="gaussian")
		tmp2 <- summary(tmp) 

		parm$clust$gamma.v <- parm$clust$beta.v <- array(0, parm$clust$G)
		parm$clust$gamma.v[selected.cols] <- 1
		parm$clust$beta.v[selected.cols] <- as.numeric(tmp2$table[(1+(1:init.rho)),1])
		parm$clust$beta0 <- as.numeric(tmp2$table[1,1])

		###############

		
		}

	#########################################

	parm$clust$beta.v <- parm$clust$beta.v[parm$clust$gamma.v==1]

	parm$clust$rho <- sum(parm$clust$gamma.v)

	#######################

	parm$missing.indx <- data$missing.indx
	parm$non.missing.indx <- data$non.missing.indx

	###############################
	######	# reset num.knots to 0 for now	
	###############################

	parm$num.knots <- 0
	# initialize to either not predictor or (predictor and 2nd order spline)
	parm$clust$gamma_3.v <- rep(0, length(parm$clust$gamma.v))
	parm$clust$rho_3 <- sum(parm$clust$gamma_3.v)

	parm <- fn.get.X.mt(parm)

	############### in case needed to set prior support for parm$sigma

	parm <- fn_varsel.assign.priors(parm, data)
	parm$sigma <- parm$prior$sigma$max #runif(n=1,parm$prior$sigma$min,parm$prior$sigma$max)
	parm$sigma.beta <- parm$prior$sigma.beta$max #runif(n=1,parm$prior$sigma.beta$min,parm$prior$sigma.beta$max)

	# gets full length (n2) parm$Z using parm$sigma asigned above
	parm <- init_fn.gen.Z(data, parm)
	
	###############################
	###############################
	###############################

	# restore num.knots to correct value from now on
	parm$num.knots <- num.knots

	# randomly pick 1 predictor and change it from 
	# 1st order to 3rd order spline -- more for debugging
	# than anything else

	pick.j <- sample(which(parm$clust$gamma.v==1), size=1)
	parm$clust$gamma_3.v[pick.j] <- 1
	parm$clust$gamma.v[pick.j] <- 0

	parm$clust$rho_3 <- sum(parm$clust$gamma_3.v)
	parm$clust$rho <- sum(parm$clust$gamma.v)

	##############################
	###### Get cluster-specific knots (fixed)
	##############################

	parm <- fn.get.knots(parm)
	
	parm <- fn.get.X.mt(parm)

	parm <- fn.gibbs_predictors_gamma(parm)
	parm <- init_fn.beta.given.gamma.sigma(data, parm)

	parm <- fn.gen.omega(data, parm)

	parm <- init_fn.gen.sigma(data, parm)

	parm$sigma.beta <- parm$sigma

	parm$test$Z <- parm$Z[data$missing.indx]
	parm$clust$test.beta.v <- parm$clust$beta.v
	parm$clust$test.beta0 <- parm$clust$beta0	

	################

	parm
	}


fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
	{if (missing(inv.sigma))
		{inv.sigma <- solve(sigma)
		}

	logdet <- as.numeric(determinant(inv.sigma, log=TRUE)$mod)
	r <- length(x)
	Q <- colSums(inv.sigma * (x-mean))
	Q <- sum(Q * (x-mean))

	val <- -r/2*log(2*pi) + logdet/2 - Q/2
	if (!log)
		{val <- exp(val)
		}

	val
	}
 


fn_varsel.quality.check <- function(parm)
	{err <- 0
	
	
  if (TRUE)
	{

	if (length(parm$clust$gamma.v) != parm$clust$G)
		{err <- 7
		}

	if (parm$clust$rho > parm$clust$G)
		{err <- 12
		}

	if (ncol(parm$clust$X.mt) != (1+parm$clust$G))
		{err <- 13
		}

	if (sum((parm$Z[data$non.missing.indx])[data$delta==0] >= data$Y[data$delta==0]) != sum(1-data$delta))
		{err <- 14
		}

	if (sum((parm$Z[data$non.missing.indx])[data$delta==1] == data$Y[data$delta==1]) != sum(data$delta))
		{err <- 15
		}

	if ((sum(is.na(parm$clust$gamma.v))>0) | (sum(is.na(parm$clust$beta.v))>0))
		{err <- 16
		}

	}
		
	err
	}




fn.init.pick.predictors <- function(parm)
	{

	# for backward compatibility
	parm$num.knots <- 0

	# CAUTION: parm$clust$X.mt has a different definition & size in main iterations
	parm$clust$X.mt <- array(,c(n2,(1+(parm$num.knots+1)*parm$clust$G)))
	parm$clust$X.mt[,1] <- 1

	parm$small.indx <- NULL
	
	for (j in 1:parm$clust$G)
		{indx.j <- which(parm$clust$c.v==j)

		knots.j <- NULL
		if (parm$num.knots>0)
			{knots.j <- as.numeric(quantile(parm$clust$A.mt[,j], prob=(1:parm$num.knots/(parm$num.knots+1))))
			}
		knots.j <- c(-Inf, knots.j, Inf)

		if (parm$clust$C.m.vec[j] > 1)
		{
    X.j.mt <- parm$all.X[,indx.j]
		a.j.v <- parm$clust$A.mt[,j]
				
		tmp.mt <- matrix(-(X.j.mt - a.j.v)^2/2, ncol=parm$clust$C.m.vec[j])
			
			L2.v <-  colSums(tmp.mt)
			
			log.prior.v <- 0 # all columns in a cluster have equal prior prob

			log.post.v <- L2.v + log.prior.v 
		
			log.post.v <- log.post.v - max(log.post.v)
			post.v <- exp(log.post.v)/ sum(exp(log.post.v))
			pick.j <- sample(indx.j, size=1, prob=post.v)
			}

		if (parm$clust$C.m.vec[j] == 1)
			{
			pick.j <- indx.j
			}

		parm$small.indx <- c(parm$small.indx, pick.j)

		expanded.X.j.mt <- NULL
		x.j.v <- parm$all.X[,pick.j]

		for (tt in 1:(parm$num.knots+1))
				{
				x.jt.v <- x.j.v*(x.j.v > knots.j[tt]) * (x.j.v <= knots.j[tt+1])
				expanded.X.j.mt <- cbind(expanded.X.j.mt, x.jt.v)
				}

		# indexing relative to parm$clust$X.mt
		start.j <- 1 + (parm$num.knots+1)*(j-1) + 1
		end.j <- 1 + (parm$num.knots+1)*j 

		parm$clust$X.mt[,start.j:end.j] <- expanded.X.j.mt

		} # END big loop "for (j in 1:parm$clust$G)"

	
	parm

	}


###################################################


fn_varsel.assign.priors <- function(parm, data)
	{
		
		parm$prior$sigma <- NULL
		parm$prior$sigma$alpha <-  1
		parm$prior$sigma$beta <-  1
		#
		
		r_sq.min <- .5
		r_sq.max <- .95
		parm$prior$sigma.sq$max <- (1-r_sq.min)*var(data$Y[data$delta==1])  
		parm$prior$sigma.sq$min <-  (1-r_sq.max)*var(data$Y[data$delta==1])  

		parm$prior$sigma$max <- sqrt(parm$prior$sigma.sq$max)
		parm$prior$sigma$min <- sqrt(parm$prior$sigma.sq$min)
		
		# 
		parm$prior$sigma.beta <- NULL
		parm$prior$sigma.beta$alpha <-  1
		parm$prior$sigma.beta$beta <-  1
		#
		######################

		sigma.beta_prop.v <- sort(1/sqrt(rgamma(n=1000,shape=parm$prior$sigma.beta$alpha,rate=parm$prior$sigma.beta$beta)))
		theta_sq.v <- (1/parm$prior$sigma.sq$min^2 + 1/sigma.beta_prop.v^2)*parm$prior$sigma.sq$min^4

		tmp <- fn.H(parm)
		H <- parm$H <- tmp[[1]]
		H_tr <- H[data$non.missing.indx, data$non.missing.indx]
		H_tr.1 <- H_tr[(data$delta==1),(data$delta==1)]
		Z.v <- data$Y[data$delta==1]

		L.v <- -length(data$non.missing.indx)/2*log(2*pi*parm$prior$sigma$min^2) -sum(Mod(eigen(H_tr.1)$val))/2*log(1+sigma.beta_prop.v^2/parm$prior$sigma$min^2)-.5/parm$prior$sigma$min^2*sum(Z.v^2)+.5*sum((H_tr.1 %*% Z.v)^2)/theta_sq.v
		parm$prior$sigma.beta$max <- sigma.beta_prop.v[which.max(L.v)]
		parm$prior$sigma.beta$min <- 0
		

		parm$prior$gamma.g$max.num.pred <- 25
		parm$prior$gamma.g$min.num.pred <- 5
		parm$prior$gamma.g$max.omega <- parm$prior$gamma.g$max.num.pred/parm$clust$G
		parm$prior$gamma.g$min.omega <- parm$prior$gamma.g$min.num.pred/parm$clust$G

		parm
	}




##################################################



init_fn.gen.Z <- function(data, parm)
	{
	parm$Z <- rep(NA,n2)

	parm <- fn.get.X.mt(parm)

	mean.v <- parm$clust$beta0 + as.vector(parm$X.mt[,-1] %*%  parm$clust$beta.v)	
	if (length(parm$missing.indx) > 0)
		{parm$Z[parm$missing.indx] <- rnorm(n=length(parm$missing.indx),mean=mean.v[parm$missing.indx], sd=parm$sigma)
		}

	# now "update" non-missing responses that are not censored
	indx1 <- intersect(parm$non.missing.indx, which(data$true$delta==1))
	parm$Z[indx1] <- data$true$Y[indx1]

	# now update non-missing responses that are censored
	indx00 <- intersect(parm$non.missing.indx, which(data$true$delta==0))
	prob.v <- pnorm(data$true$Y[indx00], mean=mean.v[indx00], sd=parm$sigma)
	u.v <- runif(n=length(indx00), min=prob.v, max=1)
	flag.v <- (prob.v > (1- 1e-5))
	#
	if (sum(flag.v) > 0)
		{parm$Z[indx00[flag.v]] <- data$true$Y[indx00[flag.v]]
		}
	parm$Z[indx00[!flag.v]] <- qnorm(u.v[!flag.v], mean=mean.v[indx00[!flag.v]], sd=parm$sigma)
		 

	parm

	}


init_fn.gen.sigma  <- function(data, parm)
	{
	###################
	# update sigma   
	###################

	a.0.v <- rep(1,n0)
	# indices rel to parm$X.mt[,-1]
	mean.v <- parm$clust$beta0*a.0.v + as.vector(matrix(parm$X.mt[parm$non.missing.indx,-1], nrow=n0) %*%  parm$clust$beta.v)
	
	resid <- parm$Z[parm$non.missing.indx] - mean.v
	shape <- parm$prior$sigma$alpha + n0/2
	rate <- parm$prior$sigma$beta + sum(resid^2)/2

	u.min <- pgamma(1/parm$prior$sigma$max^2,shape=shape, rate=rate)
	u.max <- pgamma(1/parm$prior$sigma$min^2,shape=shape, rate=rate)
	gen.u <- runif(n=1, min=u.min, max=u.max)
      
	parm$sigma <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))

	parm$maxed.sigma.flag <- FALSE

	# overwrite to avoid zeros and Inf
	if (round(u.min, dig=2) == 1) # really close to 1
		{parm$sigma <- runif(n=1, min=parm$prior$sigma$min, max=parm$prior$sigma$max)
		parm$maxed.sigma.flag <- TRUE   
		}
	if (round(u.max, dig=2) == 0) # really close to 0
		{parm$sigma <- runif(n=1, min=parm$prior$sigma$min, max=parm$prior$sigma$max) 
		parm$maxed.sigma.flag <- TRUE   
		}

	parm
	}


init_fn.beta.given.gamma.sigma <- function(data, parm)
	{

	# parm$clust$beta.v excludes the intercept
	
	
	# Call fn.get.X.mt to update parm$X.mt (which depends on gamma.v)
	# 

	parm <- fn.get.X.mt(parm)

	# always excluding intercept
	parm$clust$beta.v <- array(,(ncol(parm$X.mt)-1))

	# subset only non-missing cases
	nm.X.mt <- matrix(parm$X.mt[parm$non.missing.indx,], nrow=length(parm$non.missing.indx))
	tXX.mt <- t(nm.X.mt) %*% nm.X.mt

	# if too many covariates selected, randmly drop 
	# some to get non-singular matrices

	drop1.col.indx <- drop2.col.indx <- NULL

	# using length(parm$non.missing.indx) instead of n0
	# to protect against operator overload
	# Also when called from fn.expanded.test, n=100

	if  (ncol(tXX.mt) > length(parm$non.missing.indx))
		{# indices of alpha1.v, not parm$clust$beta.v 
		# drop 2 more than necessary just to be sure
		drop1.col.indx <- sample(1:ncol(tXX.mt), size= (ncol(tXX.mt) - length(parm$non.missing.indx) + 2) ) 
		}

	
	## sometimes, a knot is too small or too big for the covariates
	## resulting in a 0 column in the splines. Drop these columns
	## assigning beta=0 for these column coeffs (although, any 
	## finite value would work)

	drop2.col.indx <- which(colSums(abs(parm$X.mt[parm$non.missing.indx,])) < 1e-5)

	drop.col.indx <- sort(unique(c(drop1.col.indx,drop2.col.indx)))

	# making this test on basis of length rather than is.null()
	# because drop2.col.indx is not NULL even if no columns satisfy 
	# condition, but has 0 length. Same for "flag" computed below

	if  (length(drop.col.indx)>0)
		{tXX.mt <- tXX.mt[-drop.col.indx,-drop.col.indx]
		nm.X.mt <- parm$X.mt[parm$non.missing.indx,-drop.col.indx]
		}


	inv.tXX.mt <- solve(tXX.mt)

	post.var.mt <- 1/(1/parm$sigma.beta^2 + 1/parm$sigma^2) * inv.tXX.mt
	post.mean.v <- 1/parm$sigma^2/(1/parm$sigma.beta^2 + 1/parm$sigma^2) * colSums(inv.tXX.mt * colSums(nm.X.mt * parm$Z[parm$non.missing.indx]))
	
	# store for later use
	parm$beta.train <- NULL
	parm$beta.train$post.mean.v <- post.mean.v
	parm$beta.train$post.var.mt <- post.var.mt

	# note alpha1.v may contain intercept 
	# (if it isn't DROPPED in drop.col.indx)
	# but never parm$clust$beta.v

	alpha1.v <- as.vector(rmvnorm(n=1, mean=post.mean.v, sigma=post.var.mt))
	#
	# first take care of intercept

	flag <- length(drop.col.indx) > 0

	if ((!(1 %in% drop.col.indx)) & flag)
		{parm$clust$beta0 <- alpha1.v[1]
		parm$clust$beta.v[-(drop.col.indx-1)] <- alpha1.v[-1]
		parm$clust$beta.v[(drop.col.indx-1)] <- 0
		}
	if ((1 %in% drop.col.indx)  & flag)
		{parm$clust$beta0 <- 0
		drop.col.indx <- drop.col.indx[-1]
		parm$clust$beta.v[-(drop.col.indx-1)] <- alpha1.v
		parm$clust$beta.v[(drop.col.indx-1)] <- 0
		}

	if (!flag)
		{parm$clust$beta0 <- alpha1.v[1]
		parm$clust$beta.v <- alpha1.v[-1]
		}

		
	parm

	}


########################

fn.H <- function(parm)  ## calculating projection matrix H=X(X'X)^(-1)X'   after variable seletion
	{

	# Call fn.get.X.mt to update parm$X.mt (which depends on gamma.v)
	# 

	parm <- fn.get.X.mt(parm)

	# subset only non-missing cases
	nm.X.mt <- matrix(parm$X.mt[parm$non.missing.indx,], nrow=length(parm$non.missing.indx))
	tXX.mt <- t(nm.X.mt) %*% nm.X.mt

	# if too many covariates selected, randmly drop 
	# some to get non-singular matrices

	drop1.col.indx <- drop2.col.indx <- NULL

	# using length(parm$non.missing.indx) instead of n0
	# to protect against operator overload
	# Also when called from fn.expanded.test, n=100

	if  (ncol(tXX.mt) > length(parm$non.missing.indx))
		{# indices of alpha1.v, not parm$clust$beta.v 
		# drop 2 more than necessary just to be sure
		drop1.col.indx <- sample(1:ncol(tXX.mt), size= (ncol(tXX.mt) - length(parm$non.missing.indx) + 2) ) 
		}

	
	## sometimes, a knot is too small or too big for the covariates
	## resulting in a 0 column in the splines. Drop these columns
	## assigning beta=0 for these column coeffs (although, any 
	## finite value would work)

	drop2.col.indx <- which(colSums(abs(matrix(parm$X.mt[parm$non.missing.indx,],nrow=length(parm$non.missing.indx)))) < 1e-5)

	drop.col.indx <- sort(unique(c(drop1.col.indx,drop2.col.indx)))

	# making this test on basis of length rather than is.null()
	# because drop2.col.indx is not NULL even if no columns satisfy 
	# condition, but has 0 length. Same for "flag" computed below

	X.mt <- parm$X.mt

	if  (length(drop.col.indx)>0)
		{tXX.mt <- tXX.mt[-drop.col.indx,-drop.col.indx]
		X.mt <- parm$X.mt[,-drop.col.indx]
		}


	inv.tXX.mt <- solve(tXX.mt)

	# caution: unlike var selection part, this is n2 by n2,
	# not n0 by n0
	H <- X.mt %*% inv.tXX.mt %*% t(X.mt)

	list(H, drop.col.indx, inv.tXX.mt)
	}


fn.gen.Z <- function(data, parm)
	{ 
	
	parm$maxed.Z.flag <- rep(NA,length(parm$non.missing.indx))
	parm$mean_i.v <- parm$sd_i.v <- parm$mean.v <- array(,n2)

	parm <- fn.get.X.mt(parm)
	tmp.X.mt <- parm$X.mt

	#########
	# update non-missing Z's
	#########

	theta_sq <- (1/parm$sigma^2 + 1/parm$sigma.beta^2)*parm$sigma^4
	
	# H is for all n2 responses
	tmp <- fn.H(parm)
	H <- parm$H <- tmp[[1]]
	drop.col.indx <- tmp[[2]]
	inv.tXX.mt <- tmp[[3]]
	if (length(drop.col.indx) > 0)
		{tmp.X.mt <- parm$X.mt[,-drop.col.indx]
		inv.tXX.mt <- inv.tXX.mt[-drop.col.indx,-drop.col.indx]
		}
	H_11 <- H[parm$non.missing.indx,parm$non.missing.indx]
	Z.v <- parm$Z[parm$non.missing.indx]


	##########################
	# only for diagnostics
	##########################
	
	check <- FALSE
	if (check)
		{Sigma <- parm$sigma^2*diag(length(parm$non.missing.indx)) + parm$sigma.beta^2*H_11

		## solve(Sigma) equals 1/parm$sigma^2*diag(length(parm$non.missing.indx))-H_11/theta_sq)
		## determinant(Sigma,log=TRUE)$mod equals 2*length(parm$non.missing.indx)*log(parm$sigma) + sum(Mod(eigen(H_11)$val))*log(1+parm$sigma.beta^2/parm$sigma^2)	
		## dmvnorm(Z.v,mean=rep(0,length(parm$non.missing.indx)),Sigma,log=TRUE) equals -length(data$non.missing.indx)/2*log(2*pi*parm$sigma^2) -sum(Mod(eigen(H_11)$val))/2*log(1+parm$sigma.beta^2/parm$sigma^2)-.5/parm$sigma^2*sum(Z.v^2)+.5*sum((H_11 %*% Z.v)^2)/theta_sq
		}


	for (i in 1:length(parm$non.missing.indx))
		{
		## index "ccc" in full 1:n2
		ccc <- parm$non.missing.indx[i]
		delta_i <- data$true$delta[ccc]
		Y_i <- data$true$Y[ccc]

		if (delta_i==0)
			{
			h_ii <- H_11[i,i]
			H_i.v <- H_11[i,-i]

			parm$mean_i.v[ccc] <- mean.i <- sum(H_i.v * Z.v[-i])/theta_sq/(1/parm$sigma^2 - h_ii/theta_sq)
			parm$sd_i.v[ccc] <- sd.i <- 1/sqrt(1/parm$sigma^2 - h_ii/theta_sq)

			if (check)
				{mean_2.i <- as.vector(Sigma[i,-i]%*%solve(Sigma[-i,-i])%*%Z.v[-i])
				sd_2.i <- sqrt(Sigma[i,i] - as.vector(Sigma[i,-i]%*%solve(Sigma[-i,-i])%*%Sigma[i,-i]))

				if ((abs(mean_2.i-mean.i) > 1e-4) | (abs(sd_2.i-sd.i) > 1e-4))
					{stop("incorrect formula")
					}
				}
			
			prob.i <- pnorm(Y_i, mean.i, sd.i)
			flag.i <- (prob.i > (1- 1e-5))
			#
			if (flag.i)
				{Z.v[i] <- Y_i
				parm$mean.v[ccc] <- Y_i
				}
			if (!flag.i)
				{u <- runif(n=1, min=prob.i)
				Z.v[i] <- qnorm(u, mean.i, sd.i)

				# expected value of upper tail truncated normal
				# Wikipedia
				E_i <- mean.i + sd.i*dnorm((Y_i- mean.i)/sd.i)/(1-prob.i)
				parm$mean.v[ccc] <- E_i
				}
			
			}

		if (delta_i==1)
			{flag.i <- 0
			Z.v[i] <- Y_i
			parm$mean_i.v[ccc] <- parm$mean.v[ccc] <- Y_i
			parm$sd_i.v[ccc] <- 0
			}	

		parm$maxed.Z.flag[i] <- as.numeric(flag.i)
	
		}

	parm$Z[parm$non.missing.indx] <- Z.v

	# no intercept fit: 1st column of tmp.X.mt is all 1's
	# print(paste("tmp.X.mt",tmp.X.mt))
	tmp.lse <- lm(Z.v ~ 0 + tmp.X.mt[parm$non.missing.indx,])

	if (length(parm$missing.indx) > 0)
		{H_22 <- H[parm$missing.indx,parm$missing.indx]
		
		Sigma_22.1 <- parm$sigma^2*diag(length(parm$missing.indx)) + 1/(1/parm$sigma^2+1/parm$sigma.beta^2)*H_22

		# sometimes on weird occasions, eigen(Sigma_22.1)
		# gives "error code 1 from Lapack routine 'dsyevr'
		# for no good Googled reason, but laundering it 
		# in this way works and solves the problem!

		if (mode(try(eigen(Sigma_22.1)))=="character")
			{Sigma_22.1 <- solve(solve(Sigma_22.1))
			}

		beta.lse_tr <- as.vector(tmp.lse$coeff)
		# print(paste("beta.lse_tr",beta.lse_tr))
		if(length(beta.lse_tr)==1){hat.Z_test <- as.vector(tmp.X.mt[parm$missing.indx,] * beta.lse_tr)
		}else{
		hat.Z_test <- as.vector(tmp.X.mt[parm$missing.indx,] %*% beta.lse_tr)}
		mu_2.1 <- 1/parm$sigma^2/(1/parm$sigma^2+1/parm$sigma.beta^2)*hat.Z_test
		parm$mean.v[parm$missing.indx] <- mu_2.1
		
		parm$Z[parm$missing.indx] <- as.vector(rmvnorm(n=1,mean=mu_2.1, sigma=Sigma_22.1))
		}

	parm

	}


fn.get.knots <- function(parm)
	{
	
	parm$knots.mt <- matrix(apply(parm$all.X[parm$non.missing.indx,], 2, mean), nrow=parm$num.knots)

	parm
	}


fn.get.X.mt <- function(parm)
	{

	parm$p.X.mt <- parm$clust$X.mt[,-1] 	
	p.X.knots.gamma.mt <- p.X.gamma.mt <- p.X3.knots.gamma_3.mt <- p.X3.gamma_3.mt <- NULL
	
	if (parm$clust$rho > 0)
		{p.X.gamma.mt <- matrix(parm$p.X.mt[,(parm$clust$gamma.v==1)], nrow=n2)

		if (parm$num.knots > 0)
			{# parm$num.knots X parm$clust$rho
			knots.gamma.mt <- matrix(parm$knots.mt[,parm$small.indx[parm$clust$gamma.v==1]],ncol=parm$clust$rho) 
			}

		p.X.knots.gamma.mt <- array(,c(n2, (1+parm$num.knots)*parm$clust$rho))
		
		for (j in 1:parm$clust$rho)   ## needs to figure out num.knots and calculation of col.jk
			{
			col.j1 <- (parm$num.knots+1)*(j-1)+1
			p.X.knots.gamma.mt[,col.j1] <- p.X.gamma.mt[,j] 

			if (parm$num.knots > 0)
				{
				for (k in 2:(parm$num.knots+1))
					{col.jk <- ((parm$num.knots+1)*(j-1)+k)
			
					diff.jk <- p.X.gamma.mt[,j] - knots.gamma.mt[(k-1),j]
					p.X.knots.gamma.mt[,col.jk] <- diff.jk * (diff.jk > 0) 
					}
				}
			}

		} # end if


	if (parm$clust$rho_3 > 0)
		{
		p.X3.gamma_3.mt <- matrix(parm$p.X.mt[,(parm$clust$gamma_3.v==1)], nrow=n2)

		# calling this objects "knots" just for backward (and possibly forward) compatibility
		p.X3.knots.gamma_3.mt <- p.X3.gamma_3.mt

		} # end if


	# store just 1st order knots
	parm$X_1.mt <- p.X.knots.gamma.mt
	
	# store just 3rd order knots
	parm$X_3.mt <- p.X3.knots.gamma_3.mt

	parm$X.mt <- cbind(rep(1,n2), parm$X_1.mt, parm$X_3.mt)

	if (ncol(parm$X.mt) != (1+parm$clust$rho*(1+parm$num.knots)+parm$clust$rho_3) )
		{stop("Wrong number of columns in parm$X.mt")
		}

	################

	parm
	}



#####################################

fn.ortho.X.mt <- function(parm)   ## why do we need to orthogonalize, efficient calculation?
	{	
	###########################
	# Assuming intercept is always present
	# The P covariates are later orthogonalized wrt intercept 
	###########################

	# parm$E.mt will have orthonormal columns and has n0 rows

	# add intercept first of all
	tmp.X.mt <- array(1,c(n2,1))

	## CAUTION: extremely inportant to remove missing,
	# otherwise, orthogonalizing procedure will be incorrect
	
	# remove missing individuals
	tmp.X.mt <- tmp.X.mt[parm$non.missing.indx,]

	#####################################

	parm$E.mt <- qr.Q(qr(tmp.X.mt))

	# matt <- round(t(parm$E.mt) %*%parm$E.mt, dig=1) is diagonal matrix 
	# of dimension 1
	# It has only 0s and 1s on diagonal
	# sum(matt) == sum(diag(matt))

	parm
	}


fn.gibbs_predictors_gamma <- function(parm)
	{

	## store just in case
	old_old.parm <- parm

	# get parm$X_1.mt, and parm$X_3.mt
	parm <- fn.get.X.mt(parm)

	###########################
	# Assuming intercept is always present
	# The P covariates are orthogonalized wrt intercept 
	# parm$E.mt stays fixed with only intercept 
	###########################

	parm <- fn.ortho.X.mt(parm)

	###########################

	## CAUTION: extremely important to remove missing,
	# otherwise, orthogonalizing procedure will be incorrect
	
	const1 <- .5/parm$sigma^4/(1/parm$sigma^2 + 1/parm$sigma.beta^2)
	const2 <- log(1 + parm$sigma.beta^2/parm$sigma^2)
	
	h0 <-   log(1-parm$prior$gamma.g$omega-parm$prior$gamma.g$omega_3) 
	h1 <- log(parm$prior$gamma.g$omega)
	h3 <- log(parm$prior$gamma.g$omega_3)

	##################################################

	## CAUTION: Make sure both parm$clust$X.mt and parm$X.mt updated in loop
	# However, do not interactively update parm$X.mt within loop 
	# because you need to make room (delete) 
	# for (parm$num.knots+1) extra columns in the middle of 
	# parm$X.mt if parm$clust$gamma.v[j] changes 0 -> 1 (1 -> 0). Instead
	# compute parm$X.mt from scratch if parm$clust$gamma.v[j] changes 
	# from 0 -> 1. Same thing for parm$clust$gamma_3.v[j]
	
	for (j in 1:parm$clust$G) 
		{
		old.parm <- parm

		######################################################
		# check if there are too many predictors in model and if ...
		######################################################

		# Case 1: move to next j if neither regular nor spline predictor is chosen
		if ( (ncol(parm$X.mt) > 2/3*length(parm$non.missing.indx)) & (parm$clust$gamma.v[j] == 0) & (parm$clust$gamma_3.v[j] == 0) )
			{next
			}

		# Case 2: if regular (non-spline) is currently chosen, then allow only downgrading to no predictor
		# only no predictor or regular predictor
		if ( (ncol(parm$X.mt) > 2/3*length(parm$non.missing.indx)) & (parm$clust$gamma.v[j] == 1)) 
			{
			h0 <-   log(1-parm$prior$gamma.g$omega) 
			h1 <- log(parm$prior$gamma.g$omega) 
			h3 <- -Inf 
			}

		# Case 3: if spline is currently chosen, then allow all 3 options (either downgrade or stay where you are)
		if ( (ncol(parm$X.mt) > 2/3*length(parm$non.missing.indx)) & (parm$clust$gamma_3.v[j] == 1)) 
			{
			h0 <-   log(1-parm$prior$gamma.g$omega-parm$prior$gamma.g$omega_3) 
			h1 <- log(parm$prior$gamma.g$omega)
			h3 <- log(parm$prior$gamma.g$omega_3) 
			}

		#####################################################

		# -col.j has indices of parm$X.mt columns,
		# corresponding to P PREDICTORS *excluding jth cluster.
		# Depends on whether jth cluster is or isn't currently a predictor
		# The actual columns of parm$X.mt *excluding jth cluster 
		# are then stored in other.P_pred.X.mt
		
		if (parm$clust$gamma.v[j] == 0)
			{col.j <- 1
			}

		# if 1st order spline present
		if (parm$clust$gamma.v[j] == 1)
			{jjj <- which(which(parm$clust$gamma.v == 1) == j)
			col.j <- c(1, 1+((jjj-1)*(1+parm$num.knots) + (1:(parm$num.knots+1))) )
			}

		# if no spline present
		if (parm$clust$gamma_3.v[j] == 1)
			{jjj3 <- which(which(parm$clust$gamma_3.v == 1) == j)
			col.j <- c(col.j, (1+parm$clust$rho*(1+parm$num.knots) + jjj3 )) 
			}

		####
		# get remaining P predictors *excluding jth cluster
		other.P_pred.X.mt <- matrix(parm$X.mt[,-col.j], nrow=n2)
    # How to calculate the number of columns??
		if (ncol(other.P_pred.X.mt) != ( (1-parm$clust$gamma.v[j])*(1-parm$clust$gamma_3.v[j]) * (parm$clust$rho*(1+parm$num.knots) + parm$clust$rho_3)   +   parm$clust$gamma.v[j]*(1-parm$clust$gamma_3.v[j]) * ((parm$clust$rho-1)*(1+parm$num.knots) + parm$clust$rho_3)   +   (1-parm$clust$gamma.v[j])*parm$clust$gamma_3.v[j] * (parm$clust$rho*(1+parm$num.knots) + parm$clust$rho_3-1) ) )
			{stop("columns of other.P_pred.X.mt incorrect")
			}

		## CAUTION: extremely inportant to remove missing,
		# otherwise, orthogonalizing procedure will be incorrect

		tmp.X.mt <- cbind(parm$E.mt, other.P_pred.X.mt[parm$non.missing.indx,])

		if (ncol(tmp.X.mt) != 1 + ncol(other.P_pred.X.mt))
			{stop("columns of tmp.X.mt incorrect")
			}

	

		# Gram-Schmidt(tmp.X.mt)  How to carry out this??
		if (ncol(tmp.X.mt) <= n0)
			{other.E.mt <- qr.Q(qr(tmp.X.mt))
			}
		
		if (ncol(tmp.X.mt) > n0)
			{more <- ncol(tmp.X.mt) - n0
			other.E.mt <- qr.Q(qr(tmp.X.mt[,1:n0]))
			other.E.mt <- cbind(other.E.mt, array(0,c(n0,more)))
			}
		
		if (ncol(other.E.mt) != ncol(tmp.X.mt))
			{stop("columns of other.E.mt incorrect")
			}
			
		# just for checks
		matt <- round(t(other.E.mt) %*% other.E.mt, dig=1) 
		# It has only 0s and 1s on diagonal
		if (sum(matt) != sum(diag(matt)))
			{stop("matt not diagonal")
			}
		if (sum(!(diag(matt) %in% c(0,1))) > 0)
			{stop("diagonal(matt) not just 0,1")
			}

		#####################################
		#####################################

		indx.j <- which(parm$clust$c.v==j)
		X.j.mt <- matrix(parm$all.X[,indx.j], nrow=n2)   # columns that are from the jth cluster

		if (parm$num.knots > 0)
			{knots.j.mt <- matrix(parm$knots.mt[,indx.j], nrow=parm$num.knots)
			}

		## now expand X.j.mt using 1st order spline into X.j.knots.mt

		X.j.knots.mt <- array(,c(n2, (1+parm$num.knots)*parm$clust$C.m.vec[j]))
		
		for (t in 1:parm$clust$C.m.vec[j])
			{
			col.t1 <- (parm$num.knots+1)*(t-1)+1
			X.j.knots.mt[,col.t1] <- X.j.mt[,t] 

			if (parm$num.knots > 0)
				{for (k in 2:(parm$num.knots+1))
					{col.tk <- ((parm$num.knots+1)*(t-1)+k)

					diff.tk <- X.j.mt[,t]  - knots.j.mt[(k-1),t]
			
					X.j.knots.mt[,col.tk] <- diff.tk * (diff.tk > 0) 
					}
				}
			}

		## also "expand" X.j.mt using no spline into X3.j.knots.mt
		## (still calling it "knots" for backward compatibility)

		X3.j.knots.mt <- X.j.mt
				

		# clip X.j.knots.mt and X3.j.knots.mt to non-missing 
		# before orthogonalization
		# But first store full length version because will 
		# use it later to update parm$X.mt

		full.X.j.knots.mt <- X.j.knots.mt
		full.X3.j.knots.mt <- X3.j.knots.mt
		# clip it
		X.j.knots.mt <-  matrix(X.j.knots.mt[parm$non.missing.indx,], nrow=length(parm$non.missing.indx))
		X3.j.knots.mt <-  matrix(X3.j.knots.mt[parm$non.missing.indx,], nrow=length(parm$non.missing.indx))

		# Combine X.j.knots.mt and X3.j.knots.mt 	
		# into j.knots.mt
	
		j.knots.mt <- cbind(X.j.knots.mt, X3.j.knots.mt)

	
		# check 
		if (ncol(j.knots.mt) != (2+parm$num.knots)*parm$clust$C.m.vec[j])
			{stop("j.knots.mt has wrong number of columns")
			}


		## end of expansion of X.j.mt using knots into X.j.knots.mt
		## and of X3.j.mt using knots into X3.j.knots.mt
		################################################


		# projection of columns of j.knots.mt on other.E.mt
		# rows of coord.mt corres to orthogonal vectors in other.E.mt
		# cols of coord.mt corres to columns of j.knots.mt
		
		coord.mt <- t(other.E.mt) %*% j.knots.mt 
		mat <- other.E.mt %*% coord.mt 
		
		# residuals of projecting j.knots.mt on E.j.mt
		e.j.knots.mt <- j.knots.mt - mat

		# notice: columns of e.j.knots.mt are not mutually orthogonal,
		# not even within a covariate, 
		# but are each orthogonal to column space of other.E.mt
		# check:
		if (round(sum(abs(t(other.E.mt) %*% e.j.knots.mt)),dig=3) != 0)
			{stop("columns of e.j.knots.mt not orthogonal to other.E.mt")
			}

		# Now take back E.j.knots.mt and E3.j.knots.mt, 
		# the parts of X.j.knots.mt and X3.j.knots.mt
		# that are orthogonal to other.E.mt

		E.j.knots.mt <- e.j.knots.mt[,1:(parm$clust$C.m.vec[j]*(1+parm$num.knots))]
		E3.j.knots.mt <- matrix(e.j.knots.mt[,-(1:(parm$clust$C.m.vec[j]*(1+parm$num.knots)))], ncol=parm$clust$C.m.vec[j])

		## following is necessary when e.j.knots.mt is essentially 0's
		## because other.E.mt is already of full rank
		## in which case we simply generate from prior 
		## because likelihoods are the same w/ or w/o extra predictor

		skip.flag <- qr(other.E.mt)$rank == length(parm$non.missing.indx)


	if (skip.flag)
	    {L1.v <- rep(0, parm$clust$C.m.vec[j]*(1+parm$num.knots))
	    L3.v <- rep(0, parm$clust$C.m.vec[j])
	    }


	if (!skip.flag)
	    {
		# now normalize the residuals
		# so that colSums(E.j.knots.mt^2) and colSums(E3.j.knots.mt^2)
		# are all 1's or 0's
		# tmp.mt needed when e.j.knots.mt has 1 column only

		old.E.j.knots.mt <- E.j.knots.mt
		old.E3.j.knots.mt <- E3.j.knots.mt

		zero.flag.v <- round(colSums(E.j.knots.mt^2), dig=10) == 0
		if (sum(1-zero.flag.v) > 0)
			{tmp.mt <- matrix(E.j.knots.mt[,(!zero.flag.v)], nrow=length(parm$non.missing.indx))
			E.j.knots.mt[,(!zero.flag.v)] <- t(t(tmp.mt)/sqrt(colSums(tmp.mt^2)))
			}

		zero.flag.v <- round(colSums(E3.j.knots.mt^2), dig=10) == 0
		if (sum(1-zero.flag.v) > 0)
			{tmp.mt <- matrix(E3.j.knots.mt[,(!zero.flag.v)], nrow=length(parm$non.missing.indx))
			E3.j.knots.mt[,(!zero.flag.v)] <- t(t(tmp.mt)/sqrt(colSums(tmp.mt^2)))
			}


		# check:

		if (sum(!(round(colSums(E.j.knots.mt^2),dig=10) %in% c(0,1))) > 0)
			{stop("some columns of E.j.knots.mt do not have unit or ero length")
			}

		if (sum(!(round(colSums(E3.j.knots.mt^2),dig=10) %in% c(0,1))) > 0)
			{stop("some columns of E3.j.knots.mt do not have unit or ero length")
			}

		###############################################################
		#### mutually othrogonalize power 
		#### and spline columns within each covariate
		###############################################################

		##################
		## first E.j.knots.mt -> O.j.knots.mt
		##################

		O.j.knots.mt <- E.j.knots.mt

		if (parm$num.knots > 0)
			{
			for (t in 1:parm$clust$C.m.vec[j])
				{col.t <- (parm$num.knots+1)*(t-1)+ 1:(parm$num.knots+1)
				tmp.t.mt <- E.j.knots.mt[,col.t]
				qr.tmp.t.mt <- qr(tmp.t.mt)
				ortho.tmp.t.mt <- qr.Q(qr.tmp.t.mt)
				#
				# It does something weird if qr.tmp.t.mt
				# is not of full rank
				# Sometimes, numerically of full rank
				# even though one column is essentially
				# zeroes (i.e. qr()$rank cannot detect
				#
				flag1 <- qr.tmp.t.mt$rank < (parm$num.knots+1) 
				flag2 <- sum(colSums(abs(tmp.t.mt)) < 1e-5) > 0
				less.rank.flag <-  flag1 | flag2
				#
				zero.cols <- NULL
				if (less.rank.flag)
					{if (flag1)
						{zero.cols <- (qr.tmp.t.mt$rank+1):(parm$num.knots+1)
						}
					if (flag2)
						{zero.cols <- c(zero.cols,which(colSums(abs(tmp.t.mt)) < 1e-5))
						}
					zero.cols <- sort(unique(zero.cols))
					ortho.tmp.t.mt[,zero.cols] <- 0
					}
				O.j.knots.mt[,col.t] <- ortho.tmp.t.mt
				}
			# check for example: round(sum(O.j.knots.mt[,1]*O.j.knots.mt[,2]), dig=1) == 0
			}

		## notice: columns of O.j.knots.mt are mutually orthogonal
		# for the knots and powers within a covariate, 
		# but not necessarily across covariates,
		# But are each orthogonal to column space of other.E.mt
		# recheck: round(sum(abs(t(other.E.mt) %*% O.j.knots.mt)),dig=3) == 0

		if (round(sum(abs(t(other.E.mt) %*% O.j.knots.mt)),dig=3) != 0)
			{stop("columns of O.j.knots.mt not orthogonal to other.E.mt")
			}

		#
		# check:
		maxx.cor.v <- array(,parm$clust$C.m.vec[j])
		for (t in 1:parm$clust$C.m.vec[j])
			{col.t <- (parm$num.knots+1)*(t-1)+ 1:(parm$num.knots+1)
			# rather than correlations, we are only interested in dot products, i.e covariances
			cor.t.mt <- t(O.j.knots.mt[, col.t]) %*% O.j.knots.mt[, col.t]
			cor.t.mt <-round(cor.t.mt, dig=5)
			maxx.cor.v[t] <- max(abs(cor.t.mt[upper.tri(cor.t.mt)]))
			}

		if (round(max(maxx.cor.v), dig=10) != 0)
			{stop("some columns of O.j.knots.mt  within a covariate are not mutually orthogonal")
			}
			

		##################
		## next E3.j.knots.mt -> O3.j.knots.mt
		##################


		O3.j.knots.mt <- E3.j.knots.mt

		if (parm$num.knots > 0)
			{
			for (t in 1:parm$clust$C.m.vec[j])
				{col.t <- t
				tmp.t.mt <- matrix(E3.j.knots.mt[,col.t], ncol=1)
				qr.tmp.t.mt <- qr(tmp.t.mt)
				ortho.tmp.t.mt <- qr.Q(qr.tmp.t.mt)
				#
				# It does something weird if qr.tmp.t.mt
				# if not of full rank
				# Sometimes, numerically of full rank
				# even though one column is essentially
				# zeroes (i.e. qr()$rank cannot detect
				#
				flag1 <- qr.tmp.t.mt$rank < 1 
				flag2 <- sum(colSums(abs(tmp.t.mt)) < 1e-5) > 0
				less.rank.flag <-  flag1 | flag2
				#
				zero.cols <- NULL
				if (less.rank.flag)
					{if (flag1)
						{zero.cols <- (qr.tmp.t.mt$rank+1):1
						}
					if (flag2)
						{zero.cols <- c(zero.cols,which(colSums(abs(tmp.t.mt)) < 1e-5))
						}
					zero.cols <- sort(unique(zero.cols))
					ortho.tmp.t.mt[,zero.cols] <- 0
					}
				O3.j.knots.mt[,col.t] <- ortho.tmp.t.mt
				}

			# check for example: round(sum(O3.j.knots.mt[,1]*O3.j.knots.mt[,2]), dig=1) == 0
			}

		## notice: columns of O3.j.knots.mt are mutually orthogonal
		# for the knots and powers within a covariate, 
		# but not necessarily across covariates,
		# But are each orthogonal to column space of other.E.mt
		# recheck: round(sum(abs(t(other.E.mt) %*% O3.j.knots.mt)),dig=3) == 0

		if (round(sum(abs(t(other.E.mt) %*% O3.j.knots.mt)),dig=3) != 0)
			{stop("columns of O3.j.knots.mt not orthogonal to other.E.mt")
			}


		##### end check ############
		

		L1.v <- as.vector(t(O.j.knots.mt) %*% parm$Z[parm$non.missing.indx])^2
		L3.v <- as.vector(t(O3.j.knots.mt) %*% parm$Z[parm$non.missing.indx])^2

	} # end massive loop: if (!skip.flag)

		
		
		L1.mt <- matrix(L1.v, nrow=(parm$num.knots+1))
		L3.mt <- matrix(L3.v, nrow=1)

		L1.v <- colSums(L1.mt)
		L3.v <- colSums(L3.mt)

		L1.v <- -.5*(parm$num.knots+1)*const2 + const1*L1.v
		L3.v <- -.5*const2 + const1*L3.v
			
		L.v <- c(rep(0,parm$clust$C.m.vec[j]), L1.v, L3.v)

		log.prior.v <- rep(c(h0,h1,h3), each=parm$clust$C.m.vec[j])
			
		# Note: sum(exp(log.prior.v)) == parm$clust$C.m.vec[j]

		log.post.v <- L.v + log.prior.v 
		##
		log.post.v <- log.post.v - max(log.post.v)
		post.v <- exp(log.post.v)/ sum(exp(log.post.v))
		log.post.v <- log(post.v)

		gen.class <- sample(1:(3*parm$clust$C.m.vec[j]), size=1, prob=post.v) 
		gen.class <- gen.class - 1
		
		parm$clust$gamma.v[j] <- as.numeric(gen.class %/% parm$clust$C.m.vec[j] %in% c(1))
		parm$clust$gamma_3.v[j] <- as.numeric(gen.class %/% parm$clust$C.m.vec[j] %in% c(2))
		
		parm$clust$rho <- sum(parm$clust$gamma.v)
		parm$clust$rho_3 <- sum(parm$clust$gamma_3.v)

		pick.j <- gen.class %% parm$clust$C.m.vec[j] + 1
		parm$small.indx[j] <- indx.j[pick.j]
		x.j.v <- parm$all.X[,indx.j[pick.j]]

		###############################
		# IMPORTANT UPDATE parm$clust$X.mt and parm$X.mt
		###############################

		parm$clust$X.mt[,(1+j)] <- x.j.v

		no.update.flag <- ((parm$clust$gamma.v[j]==0)&(old.parm$clust$gamma.v[j]==0)) & ((parm$clust$gamma_3.v[j]==0)&(old.parm$clust$gamma_3.v[j]==0)) 
		if (!no.update.flag)
			{parm <- fn.get.X.mt(parm)
			}
	
		##################### 

		} # end for loop in j


  ##################################################

	# just to be sure
	store.X.mt <- parm$X.mt
	parm <- fn.get.X.mt(parm)

	if (round(max(abs(store.X.mt - parm$X.mt)),dig=5) > 0)
		{stop("parm$X.mt was incorrectly updated within loop")
		}

  ##################################################

	err <- fn_varsel.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC: err=",err))
		}

	parm
	}




fn.gen.omega  <- function(data, parm)
	{
	parm$clust$rho <- sum(parm$clust$gamma.v)
	parm$clust$rho_3 <- sum(parm$clust$gamma_3.v)

	prior.in.class <- 1
	prior.out.class <- 1

	##########################################
	# omega_0 marginally
	##########################################
	
	shape1 <- parm$clust$G-parm$clust$rho-parm$clust$rho_3 + prior.in.class
	shape2<- parm$clust$rho+parm$clust$rho_3 + prior.out.class
	upper.prob <- pbeta((1-parm$prior$gamma.g$min.omega), shape1=shape1, shape2=shape2)
	lower.prob <- pbeta((1-parm$prior$gamma.g$max.omega), shape1=shape1, shape2=shape2)
	gen.u <- runif(n=1, min=lower.prob, max=upper.prob)
	omega_0 <- qbeta(gen.u, shape1=shape1, shape2=shape2)
	 

	##########################################
	# (omega,omega_3) given omega_0
	##########################################

	shape1 <- parm$clust$rho_3 + prior.in.class
	shape2<- parm$clust$rho + prior.out.class
	
	parm$prior$gamma.g$omega_3 <- (1-omega_0) * rbeta(n=1, shape1=shape1, shape2=shape2)
	parm$prior$gamma.g$omega <- 1-omega_0 -parm$prior$gamma.g$omega_3
	
	parm
	}



fn.gen.both.sigma  <- function(data, parm)
	{
	##########################################
	# variance parameter of beta's
	# that is, parm$sigma.beta
	# Assuming intercept always present
	##########################################

	tmp <- fn.H(parm)
	H <- parm$H <- tmp[[1]]
	###drop.col.indx <- tmp[[2]]
	###tmp.X.mt <- parm$X.mt[,-drop.col.indx]
	H_11 <- H[parm$non.missing.indx,parm$non.missing.indx]
	Z.v <- parm$Z[parm$non.missing.indx]
	rank.H <- sum(Mod(eigen(H_11)$val))
	const <- .5*sum((H_11 %*% Z.v)^2)/parm$sigma^2

	#################

	sigma.beta_hi <- min(1.2*parm$sigma.beta, parm$prior$sigma.beta$max)
	sigma.beta_lo <- max(.8*parm$sigma.beta, parm$prior$sigma.beta$min)
	sigma.beta_prop.v <- sort(runif(n=500,min=sigma.beta_lo, max=sigma.beta_hi))

	theta_sq.v <- (1/parm$sigma^2 + 1/sigma.beta_prop.v^2)*parm$sigma^4

	L.v <- -length(data$non.missing.indx)/2*log(2*pi*parm$sigma^2) -sum(Mod(eigen(H_11)$val))/2*log(1+sigma.beta_prop.v^2/parm$sigma^2)-.5/parm$sigma^2*sum(Z.v^2)+.5*sum((H_11 %*% Z.v)^2)/theta_sq.v
	log.prior.v <- pgamma(1/sigma.beta_prop.v^2, shape= parm$prior$sigma.beta$alpha, rate= parm$prior$sigma.beta$beta, log=TRUE)
	log.post.v <- L.v + log.prior.v
	log.post.v <- log.post.v -max(log.post.v)
	post.v <- exp(log.post.v)
	u <- runif(n=1)
	
	M <- min(which(u < cumsum(post.v)/sum(post.v))) 
	parm$sigma.beta <- sigma.beta_prop.v[M]

	#################

	sigma_hi <- min(1.2*parm$sigma, parm$prior$sigma$max)
	sigma_lo <- max(.8*parm$sigma, parm$prior$sigma$min)
	sigma_prop.v <- sort(runif(n=500,min=sigma_lo, max=sigma_hi))

	theta_sq.v <- (1/sigma_prop.v^2 + 1/parm$sigma.beta^2)*sigma_prop.v^4

	L.v <- -length(data$non.missing.indx)/2*log(2*pi*sigma_prop.v^2) -sum(Mod(eigen(H_11)$val))/2*log(1+parm$sigma.beta^2/sigma_prop.v^2)-.5/sigma_prop.v^2*sum(Z.v^2)+.5*sum((H_11 %*% Z.v)^2)/theta_sq.v
	log.prior.v <- pgamma(1/sigma_prop.v^2, shape= parm$prior$sigma$alpha, rate= parm$prior$sigma$beta, log=TRUE)
	log.post.v <- L.v + log.prior.v
	log.post.v <- log.post.v -max(log.post.v)
	post.v <- exp(log.post.v)
	u <- runif(n=1)
	
	M <- min(which(u < cumsum(post.v)/sum(post.v))) 
	parm$sigma <- sigma_prop.v[M]

		
	##################

	parm
	}





fn_varsel.hyperparameters <- function(data, parm)
	{

	parm <- fn.gen.omega(data, parm)

	parm <- fn.gen.both.sigma(data, parm)

	parm

	}


############################################




fn.predict <- function(data, parm)
	{

	tmp.parm <- parm
	tmp.parm$missing.indx <- NULL
	tmp.parm$non.missing.indx <- 1:n2
	tmp.parm <- fn.gen.Z(data, tmp.parm)

	parm$test$mean.v <-  tmp.parm$mean.v

	parm

	}

########################################

fn_varsel.iter <- function(data, parm, do.gamma.prob)
	{

	err <- fn_varsel.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC 0: err=",err))
		}

	###########################################

	parm$change.gamma <- NA 

	if (as.logical(rbinom(n=1,size=1,prob=do.gamma.prob)))
	{
	
	###### UPDATE gamma.v 

	old.parm <- parm
	parm$change.gamma <- 0

	parm <- fn.gibbs_predictors_gamma(parm)

	if ((sum(parm$clust$gamma.v != old.parm$clust$gamma.v) + sum(parm$clust$gamma_3.v != old.parm$clust$gamma_3.v)) > 0)
		{parm$change.gamma <- 1
		}

	## print(c(parm$clust$rho, parm$clust$rho_3, parm$change.gamma))

	}

	###############################################

		
	err <- fn_varsel.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC 3b: err=",err))
		}
	

	##########################

	parm <- fn_varsel.hyperparameters(data, parm)

	parm <- fn.gen.Z(data, parm)

	### parm <- fn.poissonDP.hyperparm(data, parm, w=.01)

	##############

	err <- fn_varsel.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC 4: err=",err))
		}


	parm
	}



fn_varsel.mcmc <- function(All.Stuff, true, data, n.burn, n.reps, num.knots, do.gamma.prob)
	{
	# initialize

	parm <- fn_varsel.init(true, data, All.Stuff, num.knots)
	print(paste("FINISHED INITIALIZATION == ", date(),"***********"))

	init.parm <- parm
	
	err <- fn_varsel.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC at fn_varsel.init: err=",err))
		}

	for (cc in 1:n.burn)
		{old.parm <- parm

		parm <- fn_varsel.iter(data, parm, do.gamma.prob)
		
		if (TRUE)
			{# update beta's for test data
			parm <- fn.predict(data, parm)
			}

		if (cc %% n.burn == 0)
			{print(paste("BURN = ",cc,date(),"***********"))
			}
		}


	All.Stuff.Y <- NULL
	#

	#### store post-burnin entry state

	All.Stuff.Y$burnin$parm <- parm

	##############

	All.Stuff.Y$beta0.v <- array(,n.reps)
	All.Stuff.Y$G.v <- array(,n.reps)
	All.Stuff.Y$Z.v <- array(0,n2)
	# fitted line for all n2 cases
	All.Stuff.Y$mean.v <- array(0,n2)
	########################
	# for tracking Z generation-related quantities 
	# for training set individuals that are censored
	All.Stuff.Y$mean_i.v <- All.Stuff.Y$sd_i.v <- array(0,n2)
	All.Stuff.Y$sd.mean_i.v <- All.Stuff.Y$sd.sd_i.v <- array(0,n2)
	All.Stuff.Y$max.abs.Z.v <- array(0,n2)
	#########################
	All.Stuff.Y$test.mean.v <- array(0,n2) 
	All.Stuff.Y$sd.v <- array(0,n2)
	All.Stuff.Y$pred.mean.v <- array(0,n1)
	All.Stuff.Y$pred.sd.v <- array(0,n1)
	#
	All.Stuff.Y$omega.v <- array(,n.reps)
	All.Stuff.Y$rho.v <- array(,n.reps)
	All.Stuff.Y$omega_3.v <- array(,n.reps)
	All.Stuff.Y$rho_3.v <- array(,n.reps)
	All.Stuff.Y$test.rho.v <- array(,n.reps)
	All.Stuff.Y$sigma.v <- array(,n.reps)
	All.Stuff.Y$sigma.beta.v <- array(,n.reps)

	All.Stuff.Y$gamma_13.mt <- array(0,c(2,2))
	dimnames(All.Stuff.Y$gamma_13.mt) <- list(c("g1=0","g1=1"),c("g3=0","g3=1"))
	#

	# tentative record of which predictors are included
	All.Stuff.Y$binary.pred.v <- array(0,p)

	# All.Stuff.Y$d.v <- NULL

	All.Stuff.Y$tau_0.v <- All.Stuff.Y$tau.v <- All.Stuff.Y$tau_int.v <- All.Stuff.Y$G.v <- All.Stuff.Y$K.v <- array(,n.reps)
	
	# has new gamma been proposed?
	All.Stuff.Y$change.gamma.v <- array(,n.reps)

	All.Stuff.Y$maxed.Z.flag.v <- array(,n.reps)

	#
	All.Stuff.Y$fdr.v <- All.Stuff.Y$fnr.v <- All.Stuff.Y$num.pred.v <- array(,n.reps)

	All.Stuff.Y$is.pred.v <- array(0, parm$clust$G)

	for (cc in 1:n.reps)
		{parm <- fn_varsel.iter(data, parm, do.gamma.prob)

		#
		All.Stuff.Y$omega.v[cc] <- parm$prior$gamma.g$omega
		All.Stuff.Y$rho.v[cc] <- sum(parm$clust$gamma.v)
		All.Stuff.Y$omega_3.v[cc] <- parm$prior$gamma.g$omega_3
		All.Stuff.Y$rho_3.v[cc] <- sum(parm$clust$gamma_3.v)

		for (x1 in 0:1)
		for (x2 in 0:1)
			{All.Stuff.Y$gamma_13.mt[(x1+1),(x2+1)] <- All.Stuff.Y$gamma_13.mt[(x1+1),(x2+1)] + mean((parm$clust$gamma.v==x1)*(parm$clust$gamma_3.v==x2))
			}

		All.Stuff.Y$is.pred.v <- All.Stuff.Y$is.pred.v + (parm$clust$gamma.v + parm$clust$gamma_3.v > 0)

		All.Stuff.Y$sigma.v[cc] <- parm$sigma
		All.Stuff.Y$sigma.beta.v[cc] <- parm$sigma.beta

		All.Stuff.Y$maxed.Z.flag.v[cc] <- mean(parm$maxed.Z.flag)

		All.Stuff.Y$change.gamma.v[cc] <- parm$change.gamma

		All.Stuff.Y$binary.pred.v[parm$small.indx] <- All.Stuff.Y$binary.pred.v[parm$small.indx] + 1
		
		if (TRUE)
			{# update beta's for test data
			parm <- fn.predict(data, parm)
			All.Stuff.Y$test.mean.v <- All.Stuff.Y$test.mean.v + parm$test$mean.v
			All.Stuff.Y$mean.v <-  All.Stuff.Y$mean.v + parm$mean.v
			All.Stuff.Y$sd.v <-  All.Stuff.Y$sd.v + parm$mean.v^2
			All.Stuff.Y$Z.v <- All.Stuff.Y$Z.v + parm$Z
			#
			parm$pred.mean.v <- All.Stuff.Y$mean.v[data$missing.indx]
			All.Stuff.Y$pred.mean.v <- All.Stuff.Y$pred.mean.v + parm$pred.mean.v

			# for tracking Z generation-related quantities 
			# for training set individuals that are censored
			All.Stuff.Y$mean_i.v <- All.Stuff.Y$mean_i.v + parm$mean_i.v
			All.Stuff.Y$sd.mean_i.v <- All.Stuff.Y$sd.mean_i.v + parm$mean_i.v^2
			All.Stuff.Y$sd_i.v <- All.Stuff.Y$sd_i.v + parm$sd_i.v 
			All.Stuff.Y$sd.sd_i.v <- All.Stuff.Y$sd.sd_i.v + parm$sd_i.v^2 

			All.Stuff.Y$max.abs.Z.v[cc] <- max(abs(parm$Z))
			}

		if (FALSE) #(sum(abs(parm$Z)>20)>0)
			{store.parm <- parm
			stop("Large Z")
			}
		

		##########
		# For FDR, FNR calculation
		# which clusters discovered?
		#########
		
		indx1 <- union(which(parm$clust$gamma.v==1), which(parm$clust$gamma_3.v==1))
		parm$clust$discovered <- indx1
		parm$clust$not.discovered <- setdiff(1:parm$clust$G, parm$clust$discovered)
		

		##############
		## computing fdr
		###############

		# how many clusters discovered? 
		All.Stuff.Y$num.pred.v[cc] <- num.discovered <-  length(parm$clust$discovered)
		# which clusters falsely discovered? 
		falsely.discovered <- setdiff(parm$clust$discovered,true$true.pred.clusters.v)
		
		All.Stuff.Y$fdr.v[cc] <- length(falsely.discovered)/num.discovered

		
		##############
		## computing fnr
		###############

		# how many clusters not discovered? 
		num.not.discovered <- parm$clust$G - num.discovered
		# which clusters falsely not discovered? (not discovered but contain true predictors)
		falsely.not.discovered <- intersect(parm$clust$not.discovered,true$true.pred.clusters.v)

		All.Stuff.Y$fnr.v[cc] <- length(falsely.not.discovered)/num.not.discovered

		############

		if (cc %% n.reps == 0)
			{print(paste("REPS = ",cc,date(),"***********"))
			}

		} # end for loop in cc

	All.Stuff.Y$gamma_13.mt <- All.Stuff.Y$gamma_13.mt/n.reps

	All.Stuff.Y$is.pred.v <- All.Stuff.Y$is.pred.v/n.reps

	All.Stuff.Y$mean.v <-  All.Stuff.Y$mean.v/n.reps
	All.Stuff.Y$test.mean.v <- All.Stuff.Y$test.mean.v/n.reps
	All.Stuff.Y$sd.v <-  All.Stuff.Y$sd.v/n.reps - All.Stuff.Y$mean.v^2
	All.Stuff.Y$sd.v <-  sqrt(abs(round(All.Stuff.Y$sd.v, dig=10)))
	All.Stuff.Y$Z.v <- All.Stuff.Y$Z.v/n.reps

	# for tracking Z generation-related quantities 
	# for training set individuals that are censored
	All.Stuff.Y$mean_i.v <- All.Stuff.Y$mean_i.v/n.reps
	All.Stuff.Y$sd.mean_i.v <- All.Stuff.Y$sd.mean_i.v/n.reps - All.Stuff.Y$mean_i.v^2
	All.Stuff.Y$sd.mean_i.v <- sqrt(abs(round(All.Stuff.Y$sd.mean_i.v, dig=10)))
	All.Stuff.Y$sd_i.v <- All.Stuff.Y$sd_i.v/n.reps
	All.Stuff.Y$sd.sd_i.v <- All.Stuff.Y$sd.sd_i.v + parm$sd_i.v^2 
	All.Stuff.Y$sd.sd_i.v <- All.Stuff.Y$sd.sd_i.v/n.reps - All.Stuff.Y$sd_i.v^2
	All.Stuff.Y$sd.sd_i.v <- sqrt(abs(round(All.Stuff.Y$sd.sd_i.v, dig=10)))

	All.Stuff.Y$parm <- parm
	All.Stuff.Y$init.parm <- init.parm

	All.Stuff.Y$binary.pred.v <- All.Stuff.Y$binary.pred.v/n.reps
	
	All.Stuff.Y
	}

 