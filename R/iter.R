

fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
{
  
  # Computes multivariate normal density function
  # a little faster than dmvnorm function of R!
  
  if (missing(inv.sigma))
  {inv.sigma <- solve(sigma)
  }
  
  logdet <- as.numeric(determinant(inv.sigma, logarithm=TRUE)$mod)
  r <- length(x)
  Q <- colSums(inv.sigma * (x-mean))
  Q <- sum(Q * (x-mean))
  
  val <- -r/2*log(2*pi) + logdet/2 - Q/2
  if (!log)
  {val <- exp(val)
  }
  
  val
}



fn.quality.check <- function(parm)
{err <- 0

if (!is.null(parm$clust$col.nbhd))
{#if (sum(unlist(lapply(parm$clust$col.nbhd, length))) != length(parm$col.subset.I))
  #{err <- 1
  #}
}

# if (parm$tBB_flag)
#   {
#   if (sum(diag(parm$clust$tBB.mt) < 0) > 0)
# 	  {err <- 2
#     }
#   }

if (length(parm$clust$A.mt) != parm$clust$G)
{err <- 3
}

# if (ncol(parm$clust$B.mt) != (parm$clust$G+1))
# 	{err <- 4
# 	}

if ((sum(parm$clust$C.m.vec) + parm$clust$C.m0) != parm$p)
{err <- 5
}

# if (length(parm$clust$n.vec) != parm$clust$K)
# 	{err <- 6
# 	}

if (length(parm$clust$phi.v) != parm$clust$K)
{err <- 7
}

# if ((sum(parm$clust$n.vec) + parm$clust$n0) != parm$N)
# 	{err <- 8
# 	}

err
}

######################

fn.init.clusters <- function(parm)
{
  
  X.mt <- parm$X
  num.centers <- parm$G.new

  options(warn=0)
  tmp1 <- kmeans(t(X.mt), iter.max=1000, centers=num.centers, nstart=2)
  options(warn=2)
  #
  parm$clust$c.v <- tmp1$cluster
  parm$clust$G <- length(tmp1$size)
  # 
  # 
  parm$clust$C.m.vec <- array(,parm$clust$G)

  for (g in 1:parm$clust$G)
  {I.g <- (parm$clust$c.v==g)
  parm$clust$C.m.vec[g] <- sum(I.g)
  }
  
  # round(dim(X.mt)[1]/2)
  # unsupervised case
  options(warn=0)
  tmp2 <- kmeans(X.mt, iter.max=1000, centers=4, nstart=2)
  options(warn=2)
  
  # column cluster indicator
  parm$clust$r.v <- tmp2$cluster
  parm$clust$m <- parm$clust$R <- length(tmp2$size)
  parm$clust$R.m.vec <- array(,parm$clust$R)
  
  for (g in 1:parm$clust$R)
  {I.g <- (parm$clust$r.v==g)
  parm$clust$R.m.vec[g] <- sum(I.g)
  }
  
  ###########################
  
  # start from PDP model with d=0 (i.e DP)
  parm$d <- 0
  
  
  parm
}

# generate G_0 based on phi.v from initialization
# fn.G_0 <- function(alpha, parm) {
# 
#   alpha <- true_parm$clust$alpha_0
#   num_weights <- m*length(parm[[1]]$clust$phi.v)
#   betas <- rbeta(num_weights, 1, alpha)
#   remaining_stick_lengths <- c(1, cumprod(1 - betas))[1:(num_weights-1)]
#   weights <- remaining_stick_lengths * betas[1:(num_weights-1)]
#   # truncated weights add up to be 1
#   weights <- c(weights,1-round(sum(weights), digits = 14))
#   t.v <- c(1:num_weights)
#   prior.phi <- rnorm(n=length(weights),mean=parm[[1]]$clust$mu0, sd=2*parm[[1]]$clust$tau0)
#   for(mm in 1:s){
#     parm[[mm]]$clust$weights <- weights
#     parm[[mm]]$clust$t.v <- t.v
#     parm[[mm]]$clust$prior.phi <- prior.phi
#   }
#   parm
# 
# }
# 
fn.G_0 <- function(alpha, parm) {

  alpha <- true_parm$clust$alpha_0
  # weights <- parm$clust$n.c.vec / sum(parm$clust$n.c.vec)
  #
  # t.v <- c(1:num_weights)
  #
  # prior.phi <- parm$clust$phi.c.v

  for(mm in 1:s){
    parm[[mm]]$clust$weights <- parm[[mm]]$clust$n.c.vec / sum(parm[[mm]]$clust$n.c.vec)
    # parm[[mm]]$clust$t.v <- t.v
    parm[[mm]]$clust$prior.phi <- parm[[mm]]$clust$phi.c.v
  }
  parm

}

# this maybe wrong
fn.G_mt <- function(alpha, parm){
  alpha <- true_parm$clust$alpha_1
  num_weights <- length(parm$clust$prior.phi)
  
  betas <- rbeta(1, 1-parm$d, alpha+1*parm$d)
  for(i in 2:num_weights){
    betas <- c(betas, rbeta(1, 1-parm$d, alpha+i*parm$d))
  }
  remaining_stick_lengths <- c(1, cumprod(1 - betas))[1:(num_weights-1)]
  parm$clust$weights_mt <- remaining_stick_lengths * betas[1:(num_weights-1)] 
  # truncated weights add up to be 1
  parm$clust$weights_mt <- c(parm$clust$weights_mt,1-round(sum(parm$clust$weights_mt), digits = 14))
  # squeeze phi's and their weights parm$clust$weights_mt
  parm$clust$t.v_mt <- c(1:num_weights)
  parm$clust$mt.phi <- sample(parm$clust$prior.phi, num_weights, prob = parm$clust$weights, replace = TRUE)
  u.phi <- unique(parm$clust$mt.phi)
  weight_mt <- rep(0,length(u.phi))
  for(i in 1:length(u.phi)){
    weight_mt[i] <- sum(parm$clust$weights_mt[parm$clust$mt.phi==u.phi[i]])
  }
  parm$clust$mt.phi <- u.phi
  parm$clust$weights_mt <- weight_mt
  parm
  
}

fn.eda <- function(parm, data)
{
  
  parm <- fn.init.clusters(parm)
  # reintroduced on 6/29/12
  # parm$G.max <- min(parm$p/2, round(parm$clust$G*1.1))
  
  parm$Y <- array(,c(parm$clust$m,parm$clust$G)) # array(,c(parm$clust$m,parm$p))
  
  parm$clust$C.m.vec <- array(,parm$clust$G)
  
  for (g in 1:parm$clust$G){
    for (m in 1:parm$clust$m) {
  I.g <- (parm$clust$c.v==g)
  J.g <- (parm$clust$r.v==m)
  parm$clust$C.m.vec[g] <- m.g <- sum(I.g)
  parm$clust$R.m.vec[m] <- n.g <- sum(J.g)
  x.g.v <- parm$X[J.g,I.g]
  x.g.v <- mean(x.g.v)

  parm$Y[m,g] <- x.g.v
    }
  }
  
  parm$clust$C.m0 <- parm$p - sum(parm$clust$C.m.vec)
  
  parm$Y_r <- array(,c(parm$clust$m,parm$p))
  parm$clust$r.m.vec <- array(,parm$clust$m)
  
  for(g in 1:parm$clust$m){
    I.g <- (parm$clust$r.v==g)
    parm$clust$r.m.vec[g] <- m.g <- sum(I.g)
    x.g.v <- parm$X[I.g,]
    if (m.g >1){
      x.g.v <- colMeans(x.g.v)
    }
    parm$Y_r[g,] <- x.g.v
  }
  
  # parm$clust$C.m.vec <- parm$clust$r.m.vec
  # zero-th cluster?
  parm$clust$r.m0 <- parm$n2 - sum(parm$clust$r.m.vec)
  
  # TODO, name parameters in latex 
  # mass parameter of bottom layer, parm$clust$a0 in note
  parm$clust$M <- parm$a.R
  parm$clust$M0 <- .01*parm$clust$M
  
  # parm$clust$K <- data$K.max
  
  ########## Fix  DP hyperprior
  
  parm$clust$mu0 <- data$mu0
  # parm$clust$tau0 <- diff(range(as.vector(parm$X)))/6
  parm$clust$tau0 <- data$tau0
  
  #################################
  
  # parm$g ?
  # parm$g <- rep(1:parm$clust$m,each=parm$p)
  
  # parm$N <- parm$clust$m*parm$p
  
  parm$Y <- as.vector(parm$Y)
  parm$K.max <- round(true_parm$clust$alpha_0*log(length(parm$Y)))
  
  # if (computeMode$useR) {
  
  tmp <- tryCatch({
    iter.max <- ifelse((length(parm$Y) > 1000), 50, 1000)
    kmeans(parm$Y, iter.max = iter.max, centers = parm$K.max, nstart = 10,
           algorithm = "Hartigan-Wong" # TODO: MAcQueen works better?
    )}, error = function(e) {
      print("Kmeans did not converge ... using random assignment")
      cluster <- sample(1:parm$K.max, size = length(parm$Y), replace = TRUE)
      centers <- sapply(1:parm$K.max, FUN = function(x) {
        mean(parm$Y[which(cluster == x)])
      })
      list(cluster = cluster, centers = centers, size = parm$K.max)
    })
  
  
  parm$clust$s.c.v <- parm$clust$s.v <- tmp$cluster
  # all distinct phi's
  parm$clust$phi.c.v <- parm$clust$phi.v <- as.vector(tmp$centers)
  
  # might be useless, redefined later
  parm$clust$n.c.vec <-parm$clust$n.vec <- tmp$size
  
  # number of s equal to 0
  parm$clust$n0 <- 0
  
  # TODO find a correct way to update c.v's
  parm$clust$c.c.v <- parm$clust$c.v
  # change dimension of A.mt and s.mt to n2*G
  tmp.s.mt <- matrix(parm$clust$s.v, ncol = parm$clust$G)
  tmp.A.mt <- matrix(parm$clust$phi.v[tmp.s.mt],ncol= parm$clust$G)
  
  parm$clust$G <- parm$clust$K <- parm$clust$phi.list <- parm$clust$A.mt <- parm$clust$s.v <- parm$tmp.c.v <- parm$clust$c.v <- parm$clust$s.mt <- parm$clust$n.vec <- parm$clust$C.m.vec <- list()
  
  for(tt in 1:t){
    # figure out ((tt-1)*p+1):(tt*p) to accommodate different p's 
    parm$clust$c.v[[tt]] <- parm$clust$c.c.v[parm$tmp.p[[tt]]]
    parm$tmp.c.v[[tt]] <- unique(parm$clust$c.v[[tt]])
    parm$clust$G[[tt]] <- length(parm$tmp.c.v[[tt]])
    parm$clust$s.v[[tt]] <- as.vector(tmp.s.mt[,parm$tmp.c.v[[tt]]])
    parm$clust$A.mt[[tt]] <- tmp.A.mt[,parm$tmp.c.v[[tt]]]
    parm$clust$s.v[[tt]]<- as.numeric(factor(rank(parm$clust$s.v[[tt]])))
    parm$clust$s.mt[[tt]] <- matrix(parm$clust$s.v[[tt]], ncol = parm$clust$G[[tt]])
    parm$clust$K[[tt]] <- max(unique(parm$clust$s.v[[tt]]))
    parm$clust$phi.list[[tt]] <- as.vector(parm$clust$A.mt[[tt]])[unique(parm$clust$s.v[[tt]])]
    parm$clust$n.vec[[tt]] <- tabulate(parm$clust$s.v[[tt]])
    parm$clust$n.vec[[tt]] <- parm$clust$n.vec[[tt]][parm$clust$n.vec[[tt]]!=0]
  }
  
  
  
  sum.resid.sq <- 0
  
  # initialize hyperparameters
  for(tt in 1:t){
    for (mm in 1:parm$clust$m)
    {flag.v <- parm$clust$r.v == mm
    
    for (g in 1:parm$clust$G[[tt]])
    {
      X.g.mt <- parm$X[flag.v,g]
      a.g.v <- parm$clust$A.mt[[tt]][mm,g]
      resid.g.mt <- X.g.mt - a.g.v
      sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
    }
    }
  }
  
  parm$tau_int <- parm$tau <- sqrt(sum.resid.sq/parm$n2/parm$p)
  
  ###################################
  
  parm$tau_0 <- sqrt(1+parm$tau^2)
  
  # 1-parm$tau^2/var(as.vector(parm$X))
  
  # objects of full size (based on all n2 cases)
  # parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
  # 
  # if (parm$tBB_flag)
  # {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
  # }
  
  parm <- fn.assign.priors(parm, data)
  
  parm
  
  
}



fn.gen.clust <- function(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes)
{
  
  
  ###########################################
  # Missing X values
  ###########################################
  
  parm$X <- data$X
  parm$num.X.miss <- sum(is.na(parm$X))
  tmp <- which(is.na(parm$X), arr=TRUE)
  parm$X.missing.x <- tmp[,1]
  parm$X.missing.y <- tmp[,2]
  
  # Impute any missing X values by their column-specific means
  # + a small error term to guarantee non-tied values
  
  tmp.mean.v <- apply(parm$X, 2, median, na.rm=TRUE)
  tmp.sd.v <- apply(parm$X, 2, sd, na.rm=TRUE)
  if (parm$num.X.miss>0)
  { 	for (j in 1:parm$p)
  {indx.j <- is.na(parm$X[,j])
  if (sum(indx.j) > 0)
  {parm$X[indx.j,j] <- tmp.mean.v[j] + rnorm(n=sum(indx.j), sd=tmp.sd.v[j]/5)
  }
  }
  }
  
  ##################
  
  parm$G.new <- round(parm$p/2)
  parm <- fn.eda(parm, data)
  
  #################
  
  parm <- fn.hyperparameters.init(data, parm)
  
  # change it to 
  # parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes=1)
  # parm <- element_fn.update.phi.mst(parm, max.row.nbhd.size, row.frac.probes)
  # map phi_s to phi_{mst} in every t and m. 
  
  # parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
  # 
  # if (parm$tBB_flag)
  #   {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
  #   }
  
  parm
  
}

fn.init <- function(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, tBB_flag)
{
  
  parm <- NULL
  # parm$clust$G <- 10
  parm$tBB_flag <- tBB_flag
  
  # tmp <- sapply(1:t, function(x) dim(data$X[[x]]))
  # 
  # # store dimensions
  # parm$n2 <- tmp[1,] # TODO Check
  # parm$p <- tmp[2,]  # TODO Check
  
  # combine them together
  c_data <- data$X[[1]]
  if(t >= 2){
  for(i in 2:t){
    c_data <- cbind(c_data, data$X[[i]])
  }
  }
  data$X <- c_data
  
  
  parm$n2 <- dim(c_data)[1]
  # parm$n2 <- parm$clust$R
  parm$p <- dim(c_data)[2]
  
  parm$m <- m
  parm$r <- r
  
  tmp.p <- list()
  tmp.p[[1]] <- 1:cumsum(p[1,])[1]
  if(t >= 2){
  for(rr in 1:(t-1)){
    tmp.p[[rr+1]] <- (cumsum(p[1,])[rr]+1):cumsum(p[1,])[rr+1]
  }
  }
  parm$tmp.p <- tmp.p
  # mass parameter of elementwise(s) groups
  # stored later in parm$clust$M
  parm$a.R <- 1 # true$a.R
  
  # mass paramater of columns
  parm$b1 <- 0.001 # true$b1
  
  # mass paramater of column-intercept cluster
  parm$b0 <- 2.2 # true$b0
  ############################
  # For delta neighborhoods
  ############################
  
  parm$col.delta <- .1
  
  # delta-neighborhood threshold for elements
  parm$row.delta <- .1
  
  #########################################
  # generating the R- and C- clusters
  ########################################
  
  parm$shift <- 1e-04
  parm <- fn.gen.clust(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes)
  
  parm <- fn.assign.priors(parm, data)
  
  # parm$clust$s.mt <- lapply(1:t, function(x) t(parm$clust$s.mt[[x]]))
  parm$clust$s.v <- lapply(1:t, function(x) as.vector(parm$clust$s.mt[[x]]))
  # parm$clust$A.mt <- lapply(1:t, function(x) t(parm$clust$A.mt[[x]]))
  # parm$clust$G <- parm$clust$R
  # parm$clust$c.v <- parm$clust$r.v
  
  parm
  
}


fn.gen.missing.X <- function(data, parm)
{
  # impute missing X values
  
  X.mt <- data$X
  
  if (parm$num.X.miss > 0)
  {
    for (cc in 1:parm$num.X.miss)
    {i.cc <- parm$X.missing.x[cc]
    j.cc <- parm$X.missing.y[cc]
    c.cc <- parm$clust$c.v[j.cc]
    if (c.cc != 0)
    {mean.cc <- parm$clust$A.mt[i.cc, c.cc]
    }
    if (c.cc == 0)
    {mean.cc <- 1
    }
    X.mt[i.cc, j.cc] <- rnorm(n=1, mean=mean.cc, sd=parm$tau)
    }
  }
  parm$X <- X.mt
  
  parm
}


fn.standardize_orient.X <- function(parm)
{
  
  ####
  ## STANDARDIZE X columns to unit variance and zero mean
  #####
  # Do only for columns with NA's
  # For other columns, it's just a one-time calculation at the beginning of MCMC
  
  if (parm$num.X.miss > 0)
  {tmp.X <- matrix(parm$X[,parm$X.missing.y],col=parm$num.X.miss)
  mean.v <- colMeans(tmp.X)
  sd.v <- apply(tmp.X, 2, sd)
  parm$X[,parm$X.missing.y] <- t((t(tmp.X) - mean.v)/sd.v)
  }
  
  ####
  ## ORIENT X
  ####
  
  parm$X <- t(t(parm$X) * parm$clust$orient.v)
  
  parm
}


fn.assign.priors <- function(parm, data)
{
  
  parm$prior$tau <- NULL
  parm$prior$tau$alpha.tau <- 1e-2
  parm$prior$tau$beta.tau <- 1e-2
  
  # problem, 0.95 
  parm$prior$tau$max <- sqrt(.95)*sd(as.vector(data$X), na.rm=TRUE)
  parm$prior$tau$min <- 1e-10
  parm$prior$tau.sq$max <- parm$prior$tau$max^2
  parm$prior$tau.sq$min <- parm$prior$tau$min^2
  parm$prior$inv.tau.sq$max <- 1/parm$prior$tau.sq$min
  parm$prior$inv.tau.sq$min <- 1/parm$prior$tau.sq$max
  
  parm
}



########################################

fn.gen.tau.init  <- function(data, parm)
{
  ###################
  # update tau
  ###################
  
  # only covariates assigned to non-zero row and non-zero group clusters matter
  
  
  sum.resid.sq <- 0
  count <- 0
  
  # for(tt in 1:t){
  #   
  #   for (g in 1:parm$clust$m)
  #   {flag.v <- parm$clust$r.v == g
  #   z.g.v <- which(parm$clust$s.mt[[tt]][g,] > 0)
  #   
  #   if ((sum(z.g.v) > 0) & (sum(flag.v)>0))
  #   {X.g.mt <- parm$X[flag.v,z.g.v]
  #   a.g.v <- parm$clust$A.mt[[tt]][g,z.g.v]
  #   resid.g.mt <- X.g.mt - a.g.v
  #   sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  #   count <- count + sum(z.g.v)*sum(flag.v)
  #   }
  #   
  #   }
  # }
  
  for(tt in 1:t){
    for (mm in 1:parm$clust$m)
    {flag.v <- parm$clust$r.v == mm
    z.g.v <- parm$clust$c.v[[tt]] > 0
    
    if ((sum(z.g.v) > 0) & (sum(flag.v)>0))
    {X.g.mt <- parm$X[flag.v,z.g.v]
    for (g in 1:parm$clust$G[[tt]])
    {
      a.g.v <- parm$clust$A.mt[[tt]][mm,g]
      resid.g.mt <- X.g.mt - a.g.v
      sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
    }
    resid.g.mt <- X.g.mt - a.g.v
    sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
    count <- count + sum(z.g.v)*sum(flag.v)
    }
    
    }
  }
  
  shape <- parm$prior$tau$alpha + count/2
  rate <- parm$prior$tau$beta + sum.resid.sq/2
  
  u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
  u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  # overwrite to avoid zeros and Inf
  if (round(u.min, digits = 5) == 1) # really close to 1
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$min)
  }
  if (round(u.max, digits = 5) == 0) # really close to 0
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$max)
  }
  
  ###################
  # update tau_int
  ###################
  
  # only covariates assigned to intercept cluster matter
  
  sum.resid.sq <- 0
  flag.v <- parm$clust$r.v == 0
  count <- parm$clust$r.m0*parm$n2
  
  if (parm$clust$r.m0>0)
  {X.g.mt <- parm$X[,flag.v]
  a.g.v <- rep(1,parm$n2)
  resid.g.mt <- X.g.mt - a.g.v
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  }
  
  
  shape <- parm$prior$tau$alpha + count/2
  rate <- parm$prior$tau$beta + sum.resid.sq/2
  
  # shares same support as parm$tau
  u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
  u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  if (gen.u < 1e-5)
  {parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$max)
  }
  if ((1-gen.u) < 1e-5)
  {parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$min)
  }
  
  parm
}


fn.gen.tau_0  <- function(data, parm)
{
  ###################
  # update tau_0
  ###################

  sum.resid.sq <- 0
  count <- 0
  
  for (g in 1:parm$clust$G)
  {flag.v <- parm$clust$c.v == g
  z.g.v <- parm$clust$s.mt[,g] > 0
  
  if ((sum(1-z.g.v) > 0) & (sum(flag.v)>0))
  {X.g.mt <- parm$X[!z.g.v,flag.v]
  resid.g.mt <- X.g.mt
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  count <- count + sum(1-z.g.v)*sum(flag.v)
  }
  }
  
  shape <- 1 + count/2
  rate <- 1 + sum.resid.sq/2
  
  # minimum possible value of parm$tau_0 = 1.5 * maximum possible value of parm$tau
  # maximum possible value of parm$tau_0 = 3 * sd(as.vector(data$X))
  u.min <- pgamma(1/9 / var(as.vector(data$X),na.rm=TRUE),shape=shape, rate=rate)
  u.max <- pgamma(1/1.5^2/parm$prior$tau.sq$min,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau_0 <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  parm
}

fn.gen.tau  <- function(data, parm)
{
  ###################
  # update tau
  ###################
  
  # only covariates assigned to non-zero row and non-zero group clusters matter
  
  
  sum.resid.sq <- 0
  count <- 0
  
  for (g in 1:parm$clust$G)
  {flag.v <- parm$clust$c.v == g
  z.g.v <- parm$clust$s.mt[,g] > 0
  
  if ((sum(z.g.v) > 0) & (sum(flag.v)>0))
  {X.g.mt <- parm$XX[z.g.v,flag.v]
  a.g.v <- parm$clust$A.mt[z.g.v,g]
  resid.g.mt <- X.g.mt - a.g.v
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  count <- count + sum(z.g.v)*sum(flag.v)
  }
  
  }
  
  shape <- parm$prior$tau$alpha + count/2
  rate <- parm$prior$tau$beta + sum.resid.sq/2
  
  u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
  u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  # overwrite to avoid zeros and Inf
  if (round(u.min, digits = 5) == 1) # really close to 1
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$min)
  }
  if (round(u.max, digits = 5) == 0) # really close to 0
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$max)
  }
  
  ###################
  # update tau_int
  ###################
  
  # only covariates assigned to intercept cluster matter
  
  sum.resid.sq <- 0
  flag.v <- parm$clust$c.v == 0
  # count <- parm$clust$C.m0*parm$n2
  count <- parm$clust$C.m0
  
  if (parm$clust$C.m0>0)
  {X.g.mt <- parm$X[,flag.v]
  # a.g.v <- rep(1,parm$n2)
  a.g.v <- 1
  resid.g.mt <- X.g.mt - a.g.v
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  }
  
  
  shape <- parm$prior$tau$alpha + count/2
  rate <- parm$prior$tau$beta + sum.resid.sq/2
  
  # shares same support as parm$tau
  u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
  u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  parm$tau_int <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  if (gen.u < 1e-5)
  {parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$max)
  }
  if ((1-gen.u) < 1e-5)
  {parm$tau_int <- 1/sqrt(parm$prior$inv.tau.sq$min)
  }
  
  parm
}


fn.hyperparameters.init <- function(data, parm)
{
  
  # also updates update tau_int
  parm <- fn.gen.tau.init(data, parm)
  
  # parm <- fn.gen.tau_0.init(data, parm)
  
  parm
  
}

fn.hyperparameters <- function(data, parm)
{
  
  # also updates update tau_int
  parm <- fn.gen.tau(data, parm)
  
  parm <- fn.gen.tau_0(data, parm)
  
  parm
  
}


fn.funky <- function(s,t)
{# on log scale
  lgamma(s+t) - lgamma(s)
}

fn.d <- function(d, parm)
{
  
  # formula in review paper by Lijoi and Prunster
  log.lik <- sum(log(parm$b1 + (1:(parm$clust$G-1))*d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-d), (parm$clust$C.m.vec-1)))
  
  log.lik
}




fn.poissonDP.hyperparm <- function(data, parm, w=.01, max.d)
{
  
  ## update parm$d conditional on parm$b1
  ## 1/w must be an integer
  
  d.v <- seq(0,max.d,by=w)
  len <- length(d.v)
  d.v <- d.v[-len]
  len <- len-1
  
  
  log.lik.v <- sapply(d.v, fn.d, parm)
  # putting 1/2 prior mass on 0 and remaining spread uniformly on positive points in d.v
  log.p.v <- log(.5) + c(0,  rep(-log(len-1),(len-1)))
  
  log.post.v <- log.lik.v + log.p.v
  log.post.v <- log.post.v - max(log.post.v)
  post.v <- exp(log.post.v)
  post.v <- post.v/sum(post.v)
  # posterior mass is concentrate at 0 if parm$clust$G<50
  # plot(d.v, post.v, type="l")
  
  prop.d <- sample(d.v, size=1, prob=post.v)
  
  if (prop.d > 0)
  {prop.d <- runif(n=1, min=(prop.d-w), max=(prop.d+w))
  }
  
  if (prop.d != parm$d)
  {
    # MH ratio for independent proposals and
    # prior same for all d (which is true if 0 wp .5 and \in (0,max.d) wp .5)
    
    log.ratio <- fn.d(d=prop.d, parm) - fn.d(d=parm$d, parm)
    prob <- min(1, exp(log.ratio))
    flip <- rbinom(n=1, size=1, prob=prob)
    if (flip==1)
    {parm$d <- prop.d
    }
  }
  
  parm
  
}

fn.R.sq <- function(all_parm, parm, data){
  
  R_sqe  <- est <- rep(list(list()),t)
  ssto <- lapply(1:s, function(x) lapply(1:t, function(y) sum((parm$XX[[mm]][[tt]]-mean(parm$XX[[mm]][[tt]]))^2)))
  est <- lapply(1:s, function(x) lapply(1:t, function(y) t(parm$A[[mm]][[tt]][,all_parm$clust$c.v])))
  sse <- lapply(1:s, function(x) lapply(1:t, function(y) sum((parm$XX[[mm]][[tt]]-mean(est[[mm]][[tt]]))^2)))
  R_sqe <- lapply(1:s, function(x) lapply(1:t, function(y) 1-sse[[mm]][[tt]]/ssto[[mm]][[tt]]))

  for(tt in 1:t){
    est[[tt]] <- parm$clust$A.mt[[tt]]
    est[[tt]] <- t(est[[tt]][,all_parm$clust$c.v])
    
    ssto=sum((data$X-mean(data$X))^2)
    sse=sum((t(data$X)-est[[tt]])^2)
    R_sqe[[tt]] <- 1-sse/ssto
    
  }
  parm$R <- R_sqe
  parm$e <- est
  
  parm
}

fn.save <- function(parm,tmp_parm,mm,tt){
  
  # parm$clust$c.c.v[mm,,tt] <- tmp_parm$clust$c.c.v[mm,,tt]
  parm$clust$c.v[[tt]] <- tmp_parm$clust$c.v
  parm$clust$phi.list[[tt]] <- tmp_parm$clust$phi.v
  parm$clust$s.v[[tt]] <- tmp_parm$clust$s.v
  parm$clust$s.mt[[tt]] <- tmp_parm$clust$s.mt 
  parm$clust$G[[tt]]<- tmp_parm$clust$G
  parm$clust$K[[tt]] <- tmp_parm$clust$K
  parm$clust$A.mt[[tt]] <- tmp_parm$clust$A.mt
  parm$clust$n.vec[[tt]] <- tmp_parm$clust$n.vec
  parm$clust$C.m.vec[[tt]] <- tmp_parm$clust$C.m.vec
  parm$tau <- tmp_parm$tau
  parm$tau_0<- tmp_parm$tau_0
  parm$tau_int<- tmp_parm$tau_int
  
  parm
  
}

fn.save.all <- function(parm, all_parm, data){
  
  # s.mt, A.mt to R by G
  
  tmp.G <- list()
  tmp.G[[1]] <- 1:cumsum(unlist(parm$clust$G))[1]
  if(t >= 2){
  for(rr in 1:(t-1)){
    tmp.G[[rr+1]] <- (cumsum(unlist(parm$clust$G))[rr]+1):cumsum(unlist(parm$clust$G))[rr+1]
  }
  }
  all_parm$tmp.G <- tmp.G


  for(tt in 1:t){
    # change this to accommodate different p's
    parm$clust$s.mt[[tt]] <- t(all_parm$clust$s.mt[all_parm$tmp.G[[tt]],])
    # A.mt <- matrix(all_parm$clust$A.mt,nrow = dim(all_parm$X)[1])
    parm$clust$A.mt[[tt]] <- t(all_parm$clust$A.mt[all_parm$tmp.G[[tt]],])
    # matrix(A.mt[parm$tmp.p[[tt]],],nrow = length(parm$tmp.p[[tt]]))

    parm$clust$s.v[[tt]] <- as.numeric(factor(rank(c(parm$clust$s.mt[[tt]]))))
    parm$clust$K[[tt]] <- max(unique(parm$clust$s.v[[tt]]))

    a=all_parm$clust$phi.v
    parm$clust$phi.list[[tt]] <- a[sort(unique(c(parm$clust$s.mt[[tt]])))]
    parm$clust$s.mt[[tt]] <- matrix(parm$clust$s.v[[tt]],nrow = all_parm$clust$G)

    parm$clust$n.vec[[tt]] <- array(,parm$clust$K[[tt]])
    for (s in 1:parm$clust$K[[tt]])
    {parm$clust$n.vec[[tt]][s] <- sum(parm$clust$s.v[[tt]] == s)
    }
  }
  
  parm$clust$r.v <- all_parm$clust$r.v <- all_parm$clust$c.v
  parm$clust$R <- all_parm$clust$G
  parm$clust$R.m.vec <- all_parm$clust$C.m.vec
  # parm$tau <- all_parm$tau
  # parm$tau_0 <- all_parm$tau_0
  # parm$tau_int <- all_parm$tau_int
  
  parm
}


fn.comb <- function(parm, all_parm, data, mm){
  
  # XX is p by n
  # all_parm$X is G by n
  # all_parm$clust$A.mt/s.mt is G by R
  XX <- t(parm$X)
  all_parm$n2 <- dim(XX)[1]
  all_parm$p <- dim(XX)[2]
  all_parm$b1 <- 0.001
  
  tmp_c.v <- parm$clust$c.v
  c.tmp <- max(tmp_c.v[[1]])
  if(t >= 2){
  for (i in 2:t) {
    tmp_c.v[[i]] <- tmp_c.v[[i]]+c.tmp
    c.tmp <- max(unique(tmp_c.v[[i]]))
  }
  }
  
  all_parm$clust$CC.m.vec <- do.call("c", parm$clust$C.m.vec)
  all_parm$clust$cc.v <- do.call("c",tmp_c.v)
  all_parm$clust$GG <- sum(unlist(parm$clust$G))

  all_parm$X <- array(,c(all_parm$clust$GG,all_parm$p))
  
  for (g in 1:all_parm$clust$GG)
  {I.g <- (all_parm$clust$cc.v==g)
  m.g <- sum(I.g)
  x.g.v <- XX[I.g,]
  if (m.g > 1)
  {x.g.v <- colMeans(x.g.v)
  }
  all_parm$X[g,] <- x.g.v
  }
  
  all_parm$clust$C.m0 <- 0
  all_parm$clust$col.nbhd.k <- NULL
  
  all_parm$clust$phi.v <- unlist(parm$clust$phi.list)
  tmp_s.v <- parm$clust$s.v
  
  # all_parm$clust$K <- sum(unlist(parm$clust$K))
  # all_parm$clust$K <- length(all_parm$clust$phi.v)
  # save at the end of all_PDP.functions
  # all_parm$clust$R[[mm]] <- parm$clust$R
  # all_parm$clust$r.v[[mm]] <- parm$clust$r.v
  all_parm$clust$G <- parm$clust$R
  all_parm$clust$c.v <- parm$clust$r.v
  all_parm$clust$C.m.vec <- array(,all_parm$clust$G)
  
  for (g in 1:all_parm$clust$G)
  {I.g <- (all_parm$clust$c.v==g)
  all_parm$clust$C.m.vec[g] <- sum(I.g)
  }
  
  s.tmp <- max(tmp_s.v[[1]])
  if(t >=2){
  for (i in 2:t) {
    tmp_s.v[[i]] <- tmp_s.v[[i]]+s.tmp
    s.tmp <- max(unique(tmp_s.v[[i]]))
  }
  }
  # ?
  tmp_s.mt <- lapply(1:t, function(x) matrix(tmp_s.v[[x]],nrow=all_parm$clust$G))
  all_parm$clust$s.mt <- t(do.call("cbind", tmp_s.mt))

  all_parm$clust$A.mt <- t(do.call("cbind", parm$clust$A.mt))
  
  #all_parm$clust$s.mt <- apply(all_parm$clust$s.mt,2,function(x) rep(x,times=all_parm$clust$CC.m.vec))
  #all_parm$clust$A.mt <- apply(all_parm$clust$A.mt,2,function(x) rep(x,times=all_parm$clust$CC.m.vec))
  
  all_parm$clust$s.v <- c(all_parm$clust$s.mt) # unlist(tmp_s.v)

  all_parm$clust$K <- max(unique(all_parm$clust$s.v))
  all_parm$clust$n.vec <- array(,all_parm$clust$K)
  for (s in 1:all_parm$clust$K)
  {all_parm$clust$n.vec[s] <- sum(all_parm$clust$s.v == s)
  }
  # all_parm$tau <- parm$tau
  # all_parm$tau_0 <- parm$tau_0
  # all_parm$tau_int <- parm$tau_int
  
  all_parm
}


fn.est.pair <- function(x){
  
  n2 <- length(x)
  tmp.mat <- array(0,c(n2,n2))
  
  G <- max(x)
  
  for (jj in 1:G)
  {indx.jj <- which(x==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  tmp.mat
  
}

########################################

fn.iter <- function(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
	{
  # fn updates phi_{mst}
  parm <- fn.G_mt(alpha, parm)
  
  parm <- fast_PDP_fn.main(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)
  
  parm <- fn.element.phi.DP(data, parm, max.row.nbhd.size, row.frac.probes)
  
#	parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)

	# parm <- fn.hyperparameters(data, parm)

	flip <- rbinom(n=1, size=1, prob=.1)
	if (flip==1)
	  {parm <- fn.gen.missing.X(data, parm)
	}
	
	# parm$clust$B.mt <- cbind(rep(1,parm$n2), parm$clust$A.mt)
	# if (parm$tBB_flag)
	#   {parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt
	#   }

	err <- fn.quality.check(parm)
	if (err > 0)
		{stop(paste("failed QC: err=",err))
		}

	parm


	}


fn.postBurnin <- function(text, All.Stuff, offset, n.reps, data, parm.m,all_parm, dahl.flag, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
  {
  
  XX <- A <- C <- R <- rep(list(list()),t)
  min_Dahl_dist.s <- Dahl_dist.s <- est_c.v <- list()
  min_Dahl_dist.s.c <- Dahl_dist.s.c <- est_c.v.c <- rep(list(list()),t)
  min_Dahl_dist.s <- est_c.v <- lapply(1:s,function(x) min_Dahl_dist.s[[x]]=100^2)
  min_Dahl_dist.s.c <-est_c.v.c <-lapply(1:s, function(x) lapply(1:t,function(y) min_Dahl_dist.s.c[[x]][[y]]=100^2))

  all_parm.m <- list()
  mean.taxicab.v <- c()
  mean.taxicab.v.c <- array(,c(s,t))
  tmp_parm.m <- parm.m
  
  for (cc in 1:n.reps){
    
    parm.m <- fn.G_0(alpha,parm.m)
    
    for(mm in 1:s){
      
      parm <- parm.m[[mm]]
      
      for(tt in 1:t){
        data$X <- data1$X[[tt]][tmp[[mm]],]
        XX[[mm]][[tt]] <- data$X
        data <- fn.data.slice(data$X,parm,mm,tt)$data1
        tmp_parm.m[[mm]] <- fn.data.slice(data$X,parm,mm,tt)$tmp_parm
        tmp_parm <- tmp_parm.m[[mm]]
        tmp_parm$X <- data$X
        tmp_parm$p <- dim(tmp_parm$X)[2]
        
        tmp_parm <- fn.G_mt(alpha, tmp_parm)
        
        tmp_parm <- fast_PDP_fn.main(tmp_parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)
        
        tmp_parm <- fn.element.DP(data, tmp_parm, max.row.nbhd.size, row.frac.probes)
        
        tmp_parm <- fn.poissonDP.hyperparm(data, tmp_parm, w=.01, max.d=1)
        
        # tmp_parm$tau <- tau
        tmp_parm <- fn.hyperparameters(data, tmp_parm)
        
        parm <- fn.save(parm,tmp_parm,mm,tt)
        
        All.Stuff$d[[mm]][[tt]] <- tmp_parm$d
        
        # Dahl for columns
        tmp_parm$min_Dahl_dist <- min_Dahl_dist.s.c[[mm]][[tt]]
        tmp_parm$est_c.v <- est_c.v.c[[mm]][[tt]]  
        tmp_parm = fn.Dahl.c(tmp_parm, All.Stuff, dahl.flag, mm, tt)
        # mean.taxicab.v.c[mm,tt] <- mean(true_parm$clust$col.nbhd.matrix[[mm]][[tt]] != tmp_parm$dahlDist.mt)
        # if (!is.finite(mean.taxicab.v.c[mm,tt]))
        # {stop("Nan")}
        # 
        if (!dahl.flag)  # no least square allocation, calculate pairwise probability matrix
        {All.Stuff$meanDahlDist.mt.c[[mm]][[tt]] = All.Stuff$meanDahlDist.mt.c[[mm]][[tt]] + tmp_parm$dahlDist.mt
        }

        if (dahl.flag)
        {min_Dahl_dist.s.c[[mm]][[tt]]  <- tmp_parm$min_Dahl_dist
        Dahl_dist.s.c[[mm]][[tt]]  <- tmp_parm$Dahl_dist
        est_c.v.c[[mm]][[tt]]  <- tmp_parm$est_c.v
        }
        
      }
      
      aaa <- all_parm <- parm
      
      all_parm <- fn.comb(parm, all_parm, data)
      
      all_parm <- row_PDP_fn.main(all_parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)
      
      parm <- fn.save.all(parm, all_parm, data)
      
      A[[mm]] <- parm$clust$A.mt
      C[[mm]] <- parm$clust$c.v
      R[[mm]] <- parm$clust$r.v

      # Dahl for rows
      all_parm$min_Dahl_dist = min_Dahl_dist.s[[mm]]
      all_parm$est_c.v <- est_c.v[[mm]]
      all_parm = fn.Dahl(all_parm, All.Stuff, dahl.flag, mm)
      # mean.taxicab.v[mm] <- mean(true_parm$clust$row.nbhd.matrix[[mm]] != all_parm$dahlDist.mt)
      # if (!is.finite(mean.taxicab.v[mm]))
      # {stop("Nan")}

      if (!dahl.flag)  # no least square allocation, calculate pairwise probability matrix
      {All.Stuff$meanDahlDist.mt[[mm]] = All.Stuff$meanDahlDist.mt[[mm]] + all_parm$dahlDist.mt
      }

      if (dahl.flag)
      {min_Dahl_dist.s[[mm]]  <- all_parm$min_Dahl_dist
      Dahl_dist.s[[mm]]  <- all_parm$Dahl_dist
      est_c.v[[mm]]  <- all_parm$est_c.v
      }
      
      # store 
      parm.m[[mm]] <- parm
      all_parm.m[[mm]] <- all_parm
      
    }
    
    R_sqe  <- est <- rep(list(list()),t)
    ssto <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-mean(XX[[x]][[y]]))^2)))
    est <- lapply(1:s, function(x) lapply(1:t, function(y) A[[x]][[y]][R[[x]],C[[x]][[y]]]))
    sse <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-est[[x]][[y]])^2)))
    R_sqe <- lapply(1:s, function(x) lapply(1:t, function(y) 1-sse[[x]][[y]]/ssto[[x]][[y]]))
    All.Stuff$R_sq[[cc]] <- R_sqe

  All.Stuff$G.v[[cc+offset]] <- parm$clust$G
  All.Stuff$K.v[[cc+offset]] <- parm$clust$K
  All.Stuff$tau.v[[cc+offset]] <- parm$tau
  #All.Stuff$tau_0.v[[cc+offset]] <- parm$tau_0
  #All.Stuff$tau_int.v[[cc+offset]] <- parm$tau_int
  All.Stuff$d.v[[cc+offset]] <- All.Stuff$d
  #All.Stuff$row.flip.v[[cc+offset]]  <- parm$clust$row.flip
  #All.Stuff$col_new_clust.v[[cc+offset]]  <- parm$clust$col.new.flag
  #All.Stuff$col_flip.v[[cc+offset]]  <- parm$clust$col.mh.flip
  #All.Stuff$col_exit.v[[cc+offset]]  <- parm$clust$col.mh.exit
  # All.Stuff$mean.taxicab.v[[cc+offset]] <- mean.taxicab.v
  # All.Stuff$mean.taxicab.v.c[[cc+offset]] <- mean.taxicab.v.c
  
  
  # if (dahl.flag)
  # {All.Stuff$runningMinDahlDist.v[[cc]]  <- min_Dahl_dist.s
  # All.Stuff$runningMinDahlDist.v.c[[cc]] <- min_Dahl_dist.s.c
  # All.Stuff$dahlDist.v[[cc]]  <- Dahl_dist.s
  # All.Stuff$dahlDist.v.c[[cc]]  <- Dahl_dist.s.c
  # }
  
  if (cc %% 100 == 0)
    {print(paste(text, "REPS = ",cc,date(),"***********"))
    }
  
  } # END FOR LOOP

  if (!dahl.flag)
  {All.Stuff$meanDahlDist.mt = lapply(1:s,function(x) All.Stuff$meanDahlDist.mt[[x]]/n.reps)
  All.Stuff$meanDahlDist.mt.c = lapply(1:s, function(x) lapply(1:t, function(y) All.Stuff$meanDahlDist.mt.c[[x]][[y]]/n.reps))
  }

  # change this
  if (dahl.flag)
  {
  All.Stuff$est_c.v  <- est_c.v
  All.Stuff$min_Dahl_dist = min_Dahl_dist.s
  All.Stuff$est_c.v.c  <- est_c.v.c
  All.Stuff$min_Dahl_dist.c = min_Dahl_dist.s.c

  # tmp.mat <- lapply(1:s, function(x) fn.est.pair(All.Stuff$est_c.v[[x]]))
  # All.Stuff$mean.taxicab.v.e <- lapply(1:s,function(x) mean(true_parm$clust$row.nbhd.matrix[[x]] != tmp.mat[[x]]))
  # 
  # tmp.mat <- lapply(1:s, function(x) lapply(1:t, function(y) fn.est.pair(All.Stuff$est_c.v.c[[x]][[y]])))
  # All.Stuff$mean.taxicab.v.c.e <- lapply(1:s,function(x) lapply(1:t, function(y) mean(true_parm$clust$col.nbhd.matrix[[x]][[y]] != tmp.mat[[x]][[y]])))
  
  # R_sqe  <- est <- rep(list(list()),t)
  # ssto <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-mean(XX[[x]][[y]]))^2)))
  # All.Stuff$est.v <- est <- lapply(1:s, function(x) lapply(1:t, function(y) A[[x]][[y]][All.Stuff$est_c.v[[x]],All.Stuff$est_c.v.c[[x]][[y]]]))
  # sse <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-est[[x]][[y]])^2)))
  # R_sqe <- lapply(1:s, function(x) lapply(1:t, function(y) 1-sse[[x]][[y]]/ssto[[x]][[y]]))
  # 
  # All.Stuff$R_sq_e <- R_sqe

  }
  
  tmp = list(parm.m, All.Stuff, all_parm.m)
  names(tmp) = c("parm.m", "All.Stuff", "all_parm.m")
  
  tmp
  }




fn.mcmc <- function(true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, dahl.flag=FALSE,
                    tBB_flag=FALSE)
	{
  
  init.parm <- parm.m <- all_parm.m <- parm <- NULL
  data$mu0 <- mean(unlist(data$X))
  data$tau0 <- sd(unlist(data$X))
  tmp <- tmp_tau<- list()
  tmp[[1]] <- 1:cumsum(r)[1]
  if(s >1){
  for(rr in 1:(s-1)){
    tmp[[rr+1]] <- (cumsum(r)[rr]+1):cumsum(r)[rr+1]
  }
  }
  # initialize using un_PDP initialization, store result in parm.m[[mm]]
  for(mm in 1:s){
    data$X <- lapply(1:t, function(x) data1$X[[x]][tmp[[mm]],])
  parm.m[[mm]] <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, true_parm, tBB_flag)
  init.parm[[mm]] <- parm.m[[mm]]
  }

	text="BURNIN..."
	
	data <- NULL
	tmp_parm.m <- parm.m
	
  All.Stuff <- NULL

  ccc.v <- array(,c(s,dim(init.parm[[1]]$X)[1],n.burn))
  G <- array(,c(s,n.burn))
	XX <- A <- C <- R <- rep(list(list()),t)
	
	 for (cc in 1:n.burn){
	   
	   parm.m <- fn.G_0(alpha,parm.m)
	   
	   for(mm in 1:s){
	     
	     parm <- parm.m[[mm]]

	   for(tt in 1:t){
	     data$X <- data1$X[[tt]][tmp[[mm]],]
	     XX[[mm]][[tt]] <- data$X
	     data <- fn.data.slice(data$X,parm,mm,tt)$data1
	     tmp_parm.m[[mm]] <- fn.data.slice(data$X,parm,mm,tt)$tmp_parm
	     tmp_parm <- tmp_parm.m[[mm]]
	     tmp_parm$X <- data$X
	     tmp_parm$p <- dim(tmp_parm$X)[2]
	     
	     tmp_parm <- fn.G_mt(alpha, tmp_parm)
	     
	     tmp_parm <- fast_PDP_fn.main(tmp_parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)
	     
	     tmp_parm <- fn.element.DP(data, tmp_parm, max.row.nbhd.size, row.frac.probes)
	     
	     tmp_parm <- fn.poissonDP.hyperparm(data, tmp_parm, w=.01, max.d=1)
	     
	     # tmp_parm$tau <- tau
	     tmp_parm <- fn.hyperparameters(data, tmp_parm)
	     
	     parm <- fn.save(parm,tmp_parm,mm,tt)
	     
	     # All.Stuff$d[[mm]][[tt]] <- tmp_parm$d
	     
	   }
	 
	  aaa <- all_parm <- parm

	  all_parm <- fn.comb(parm, all_parm, data)

	  all_parm <- row_PDP_fn.main(all_parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)

	  parm <- fn.save.all(parm, all_parm, data)

	  ccc.v[mm,,cc] <- all_parm$clust$c.v
	  G[mm,cc] <- all_parm$clust$G
	  # image.plot(t(parm$clust$A.mt[[2]][parm$clust$r.v,parm$clust$c.v[[2]]]))
	  # image.plot(t(XX[[2]][[2]]))
	  A[[mm]] <- parm$clust$A.mt
	  C[[mm]] <- parm$clust$c.v
	  R[[mm]] <- parm$clust$r.v
 
	  parm.m[[mm]] <- parm
	 
	   }
	   
	   R_sqe  <- est <- rep(list(list()),t)
	   ssto <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-mean(XX[[x]][[y]]))^2)))
	   est <- lapply(1:s, function(x) lapply(1:t, function(y) A[[x]][[y]][R[[x]],C[[x]][[y]]]))
	   sse <- lapply(1:s, function(x) lapply(1:t, function(y) sum((XX[[x]][[y]]-est[[x]][[y]])^2)))
	   R_sqe <- lapply(1:s, function(x) lapply(1:t, function(y) 1-sse[[x]][[y]]/ssto[[x]][[y]]))
	   All.Stuff$R_sq[[cc]] <- R_sqe

	   tmp_tau[[cc]] <- parm$tau
	   
	  if (cc %% 100 == 0)
	  {print(paste(text, "REPS = ",cc,date(),"***********"))
	  }
	  
		}
	
	##########################################
	## DEFINE OBJECTS
	##########################################

	All.Stuff <- NULL
	# redefine these
	All.Stuff$d.v <- All.Stuff$tau_0.v <- All.Stuff$tau.v <- All.Stuff$tau_int.v <- All.Stuff$G.v <- All.Stuff$K.v <- list()
	All.Stuff$row.flip.v  <- list()
	All.Stuff$nbhd_max <- All.Stuff$col_new_clust.v  <- All.Stuff$col_exit.v <- All.Stuff$col_flip.v  <- list()

	# this is only for n.reps (runs 1 or 2)
	All.Stuff$meanDahlDist.mt <- lapply(1:s, function(x) array(0,c(r[x],r[x])))
	# for(ss in 1:s){
	# All.Stuff$meanDahlDist.mt[[ss]] <- array(0,c(r[ss],r[ss]))
	# }
	All.Stuff$meanDahlDist.mt.c <- lapply(1:s, function(x) lapply(1:t, function(y) array(0,c(p[x,y],p[x,y]))))
	
	All.Stuff$d <- rep(list(list()),s)
	All.Stuff$R_sq.v <- All.Stuff$est.v <- list()

	# this is only for n.reps (runs 1 or 2)
	All.Stuff$mean.taxicab.v <- All.Stuff$runningMinDahlDist.v  <- All.Stuff$dahlDist.v <- list()
	All.Stuff$mean.taxicab.v.c <- All.Stuff$runningMinDahlDist.v.c <- All.Stuff$dahlDist.v.c <- list()
	##########################################
	## POST-MCMC RUN 1
	##########################################
	
	text="POST--BURNIN RUN 1..."
	
	# problem when dahl.flag=FALSE
	tmp_1 = fn.postBurnin(text, All.Stuff, offset=0, n.reps, data, parm.m,all_parm, dahl.flag=FALSE, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
  parm.m = tmp_1$parm.m
  All.Stuff = tmp_1$All.Stuff
  all_parm.m = tmp_1$all_parm.m

  # parm$min_Dahl_dist <- parm$Dahl_dist <- rep(list(list()),t)
  
	###
	
	##########################################
	## POST-MCMC RUN 2: ALSO COMPUTE LEAST SQUARES ALLOCATION
	##########################################
	
	text="RUN 2: DAHL..."
	
	tmp_1 = fn.postBurnin(text, All.Stuff, offset=n.reps, n.reps, data, parm.m,all_parm, dahl.flag=TRUE, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm)
	parm.m = tmp_1$parm.m
	All.Stuff = tmp_1$All.Stuff
	all_parm.m = tmp_1$all_parm.m
	
	#######################
	
	All.Stuff$parm.m <- parm.m
	All.Stuff$all_parm.m <- all_parm.m
	# All.Stuff$init.parm <- init.parm

	All.Stuff
	}

