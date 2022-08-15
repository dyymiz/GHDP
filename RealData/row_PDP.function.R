
row_PDP_fn.compact.consistency.check <- function(parm)
{err <- 0

if (max(sort(unique(parm$clust$c.v))) != parm$clust$G)
{err <- 1
}

if (sum(parm$clust$C.m.vec==0) > 0)
{err <- 1.5
}

err <- row_PDP_fn.consistency.check(parm)


err
}

row_PDP_fn.check.nbhd <- function(parm)
{
  max_dist.v <- array(,length(parm$clust$col.nbhd.k))
  
  for (zz in 1:length(parm$clust$col.nbhd.k))
  {
    v1 <- parm$clust$col_post.prob.mt[,parm$clust$col.nbhd.k[zz]]
    m1 <- matrix(parm$clust$col_post.prob.mt[,parm$clust$col.nbhd[[zz]]],ncol=length(parm$clust$col.nbhd[[zz]]))
    max_dist.v[zz] <- max(2*(1-colSums(sqrt(v1*m1))))
  }
  
  parm$clust$nbhd_max_dist <- max(max_dist.v)
  
  parm
  
}

row_PDP_fn.consistency.check <- function(parm)
{err <- 0


if (sum(parm$clust$C.m.vec) + parm$clust$C.m0 != parm$p)
{err <- 4
}

if (length(parm$clust$n.vec) != parm$clust$K)
{err <- 9
}

if (length(parm$clust$phi.v) != parm$clust$K)
{err <- 10
}


err
}


row_PDP_fn.swap.clusters <- function(parm, g1, g2)
{
  
  ####################################################
  # swap the group labels g1 and g2
  ####################################################
  
  ind1 <- parm$clust$c.v == g1
  ind2 <- parm$clust$c.v == g2
  parm$clust$c.v[ind1] <- g2
  parm$clust$c.v[ind2] <- g1
  
  buffer <- parm$clust$s.mt[,g1]
  parm$clust$s.mt[,g1] <- parm$clust$s.mt[,g2] # HOT 3
  parm$clust$s.mt[,g2] <- buffer
  
  buffer <- parm$clust$beta.v[g1]
  parm$clust$beta.v[g1] <- parm$clust$beta.v[g2]
  parm$clust$beta.v[g2] <- buffer
  
  buffer <- parm$clust$gamma.v[g1]
  parm$clust$gamma.v[g1] <- parm$clust$gamma.v[g2]
  parm$clust$gamma.v[g2] <- buffer
  
  buffer <- parm$clust$C.m.vec[g1]
  parm$clust$C.m.vec[g1] <- parm$clust$C.m.vec[g2]
  parm$clust$C.m.vec[g2] <- buffer
  
  buffer <- parm$clust$small.indx[g1]
  parm$clust$small.indx[g1] <- parm$clust$small.indx[g2]
  parm$clust$small.indx[g2] <- buffer
  
  buffer <- parm$clust$order.v[g1]
  parm$clust$order.v[g1] <- parm$clust$order.v[g2]
  parm$clust$order.v[g2] <- buffer
  
  #####################
  
  buffer <- parm$clust$A.mt[,g1]
  parm$clust$A.mt[,g1] <- parm$clust$A.mt[,g2]
  parm$clust$A.mt[,g2] <- buffer
  
  # buffer <- parm$clust$B.mt[,(g1+1)]
  # parm$clust$B.mt[,(g1+1)] <- parm$clust$B.mt[,(g2+1)]
  # parm$clust$B.mt[,(g2+1)] <- buffer
  # 
  # if (parm$tBB_flag)
  # {
  # 
  #   # first swap columns
  #   buffer <- parm$clust$tBB.mt[,(g1+1)]
  #   parm$clust$tBB.mt[,(g1+1)] <- parm$clust$tBB.mt[,(g2+1)] # HOT 3
  #   parm$clust$tBB.mt[,(g2+1)] <- buffer
  #   # then swap rows
  #   buffer <- parm$clust$tBB.mt[(g1+1),]
  #   parm$clust$tBB.mt[(g1+1),] <- parm$clust$tBB.mt[(g2+1),]
  #   parm$clust$tBB.mt[(g2+1),] <- buffer
  # 
  # } # end if (parm$tBB_flag)
  
  parm
}

row_PDP_fn.clip.clusters <- function(parm, keep)
{
  
  parm$clust$s.mt <- parm$clust$s.mt[,keep]
  
  parm$clust$beta.v <- parm$clust$beta.v[keep]
  parm$clust$gamma.v <- parm$clust$gamma.v[keep]
  parm$clust$C.m.vec <- parm$clust$C.m.vec[keep]
  parm$clust$small.indx <- parm$clust$small.indx[keep]
  parm$clust$order.v <- parm$clust$order.v[keep]
  
  parm$clust$A.mt <- parm$clust$A.mt[,keep]
  
  indx2 <- c(1,(keep+1))
  # parm$clust$B.mt <- parm$clust$B.mt[,indx2]
  # 
  # if(parm$tBB_flag)
  # {parm$clust$tBB.mt <- parm$clust$tBB.mt[indx2,indx2]
  # }
  
  parm
}



###########################################################


row_PDP_fn.log.lik <- function(gg, x.mt, parm, colSums=TRUE)
{
  if (gg > 0)
  {a2.v <- parm$clust$A.mt[,gg]
  z.g.v <- parm$clust$s.mt[,gg] > 0
  log.lik.v <- rep(0,ncol(x.mt)) ## HOT
  
  if (sum(z.g.v) > 0)
  {a2.1.v <- parm$clust$A.mt[z.g.v,gg]
  small.X.1 <- matrix(x.mt[z.g.v,], ncol=ncol(x.mt)) # HOT
  tau <- c(parm$tau/sqrt(parm$clust$CC.m.vec))
  
  if (colSums)
  {tmp <- colSums(-.5*(small.X.1 - a2.1.v)^2/tau^2)
  # tmp <- colSums(-.5*(small.X.1 - a2.1.v)^2)
  } # end 	if (colSums)
  
  if (!colSums)
  {tmp <- sum(-.5*(small.X.1 - a2.1.v)^2/tau^2)
  # tmp <- sum(-.5*(small.X.1 - a2.1.v)^2)
  } # end if (!colSums)
  
  log.lik.v <- log.lik.v + tmp + sum(z.g.v)*(-.5*log(2*pi))-sum(log(tau))
  # log.lik.v <- log.lik.v + tmp/parm$tau^2 + sum(z.g.v)*(-.5*log(2*pi)-log(parm$tau))
  } # end if (sum(z.g.v) > 0)
  
  if (sum(1-z.g.v) > 0)
  {small.X.0 <- matrix(x.mt[!z.g.v,], ncol=ncol(x.mt))
  if (colSums)
  {tmp <- colSums(-.5*small.X.0^2)
  }
  if (!colSums)
  {tmp <- sum(-.5*small.X.0^2)
  }
  log.lik.v <- log.lik.v + tmp/parm$tau_0^2 + sum(1-z.g.v)*(-.5*log(2*pi)-log(parm$tau_0))
  } # end if (sum(1-z.g.v) > 0)
  } # end if (gg > 0)
  
  if (gg == 0)
  {a2.v <- rep(1,parm$n2)
  small.X <- x.mt
  if (colSums)
  {tmp <- colSums(-.5*(small.X - a2.v)^2)
  }
  if (!colSums)
  {tmp <- sum(-.5*(small.X - a2.v)^2)
  } # end if (gg == 0)
  
  log.lik.v <-  tmp/parm$tau_int^2 - parm$n2*.5*log(2*pi)-parm$n2*log(parm$tau_int)
  }
  
  log.lik.v
}



row_PDP_fn.nbhd <- function(relative_I, parm, max.col.nbhd.size)
{
  if (length(relative_I)>1)
  {relative_k <- sample(relative_I, size=1)
  }
  if (length(relative_I)==1)
  {relative_k <- relative_I
  }
  
  post.prob.mt <- parm$clust$col_subset_post.prob.mt
  
  tmp1.mt <- matrix(post.prob.mt[,relative_I], ncol = length(relative_I)) ## HOT
  tmp2.v <- post.prob.mt[,relative_k]
  tmp3.mt <- sqrt(tmp1.mt * tmp2.v) ## HOT
  H.v <-  2 * (1 - colSums(tmp3.mt))
  
  cutoff <- parm$col.delta
  flag.v <- which(H.v <= cutoff)
  relative_I.k <- relative_I[flag.v]
  
  if (length(relative_I.k) > max.col.nbhd.size) {
    relative_I.k <- relative_I[rank(H.v, ties="random") <= max.col.nbhd.size]
  }
  
  relative_I.k <- sort(relative_I.k)
  
  relative_I <- sort(setdiff(relative_I, relative_I.k))
  relative_I <- sort(relative_I)
  
  list(relative_k, relative_I.k, relative_I)
  
}


row_PDP_fn.post.prob.and.delta <- function(parm, max.col.nbhd.size, col.frac.probes) {
  
  # change this to 1:parm$clust$G
  col.subset <- 1:parm$p
  
  
  ################################################
  ### Compute pmf of cluster variables w_1,...,w_p
  ###############################################
  
  prior.prob.v <- c(parm$clust$C.m0, parm$clust$C.m.vec)
  small <- 1e-3 # compared to 1
  prior.prob.v[prior.prob.v < small] <- small
  
  subset_log.ss.mt <- array(, c((parm$clust$G + 1), length(col.subset)))
  
  for (gg in 0:parm$clust$G)
  {subset_log.ss.mt[(gg+1),] <- row_PDP_fn.log.lik(gg, x.mt = parm$X[,col.subset], parm)
  }
  
  subset_log.ss.mt <- subset_log.ss.mt + log(prior.prob.v)
  
  maxx.v <- apply(subset_log.ss.mt, 2, max)
  
  subset_log.ss.mt <- t(t(subset_log.ss.mt) - maxx.v)
  subset_ss.mt <- exp(subset_log.ss.mt)
  
  col.sums.v <- colSums(subset_ss.mt)
  subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)
  
  # TODO Ask GS if second normalization is really necessary
  
  # replace zeros by "small"
  small2 <- 1e-5
  subset_ss.mt[subset_ss.mt < small2] <- small2
  
  # again normalize
  col.sums.v <- colSums(subset_ss.mt)
  subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)
  
  parm$clust$col_post.prob.mt <- array(,c((parm$clust$G+1), parm$p))
  parm$clust$col_post.prob.mt[,col.subset] <- subset_ss.mt
  
  dimnames(parm$clust$col_post.prob.mt) <- list(0:parm$clust$G, 1:parm$p)
  
  parm$clust$col_subset_post.prob.mt <- subset_ss.mt
  dimnames(parm$clust$col_subset_post.prob.mt) <- list(0:parm$clust$G, 1:length(col.subset))
  
  #########################################
  ### now compute the delta-neighborhoods
  #########################################
  
  
  parm$clust$col.nbhd <- NULL
  parm$clust$col.nbhd.k <- NULL
  relative_I <- 1:length(col.subset)
  
  while (length(relative_I) >= 1) {
    tmp <- row_PDP_fn.nbhd(relative_I, parm, max.col.nbhd.size) # HOT inside function
    relative_k <- tmp[[1]]
    relative_I.k <- tmp[[2]]
    relative_I <- tmp[[3]]
    #
    parm$clust$col.nbhd <- c(parm$clust$col.nbhd, list(col.subset[relative_I.k]))
    parm$clust$col.nbhd.k <- c(parm$clust$col.nbhd.k, col.subset[relative_k])
  }
  
  ### SG: how good are these nbhds?
  
  parm <- row_PDP_fn.check.nbhd(parm)
  
  
  
  ## END
  
  parm
}


###########################################################

row_PDP_fn.gibbs <- function(k, parm, data)
{	k <- parm$k

err <- row_PDP_fn.consistency.check(parm)
if (err > 0)
{stop(paste("GIBBS - 0: failed consistency check: err=",err))
}

old.c.k <- parm$clust$c.v[k]

###############

if (old.c.k > 0)
{parm$clust$C.m.vec[old.c.k] <- parm$clust$C.m.vec[old.c.k] - 1
}
if (old.c.k == 0)
{parm$clust$C.m0 <- parm$clust$C.m0 - 1
}

x.mt <- matrix(parm$X[,k], ncol=1)

# intercept cluster or any existing cluster
L.v <- sapply(0:parm$clust$G, row_PDP_fn.log.lik, x.mt, parm) # HOT


#######################################################
### emptied clusters are gone forever under Gibbs sampling
#######################################################

emptied.indx <- which(parm$clust$C.m.vec==0)
new.G <- parm$clust$G - length(emptied.indx)

if (length(emptied.indx) >0)
{
  
  new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
  new.n.vec <- array(,parm$clust$K)
  for (pp in 1:parm$clust$K) {
    new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
  }
  new.n0 <- sum(new.s.mt == 0) # HOT 3
  
  emptied.s.indx <- which(new.n.vec == 0)
  new.K <- parm$clust$K - length(emptied.s.indx)
  
}

if (length(emptied.indx) ==0)
{
  new.s.mt <- parm$clust$s.mt
  new.n.vec <- parm$clust$n.vec
  emptied.s.indx <- which(new.n.vec==0)
  new.K <- parm$clust$K
  new.n0 <- parm$clust$n0
}

## generate auxilliary P vector

tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
tmp.alpha <- tmp.M+new.n.vec
tmp.alpha[emptied.s.indx] <- 0
P.aux <- rgamma(parm$clust$K+1,c(0,tmp.alpha),1)
# rgamma(parm$clust$K+1,c(parm$clust$M0+new.n0,tmp.alpha),1)
P.aux <- P.aux/sum(P.aux)

## marginal likelihood of new cluster

marg.log.lik.v <- array(,length(x.mt))
for (tt in 1:length(x.mt))
{
  tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) # HOT 3
  tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
  
  marg.log.lik.v[tt] <- log(sum(tmp.lik.v*P.aux))
}
marg.log.lik <- sum(marg.log.lik.v)


L.v <- c(L.v, marg.log.lik)

log.prior.v <- array(NA, (2+parm$clust$G))

# allow -Inf's in Gibbs sampling log-prior (just emptied clusters)
if (length(emptied.indx) >0)
{log.prior.v[-(emptied.indx+1)] <- log(c(0, (parm$clust$C.m.vec[-emptied.indx]-parm$d), 0)) # (parm$b1+new.G*parm$d)))
# log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), (parm$b1+new.G*parm$d)))
log.prior.v[emptied.indx+1] <- -Inf
}

if (length(emptied.indx) ==0)
{log.prior.v <- log(c(0, (parm$clust$C.m.vec-parm$d), 0)) # (parm$b1+new.G*parm$d)))
# log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
}

tmp2 <- log.prior.v + L.v
maxx <- max(tmp2)
tmp2 <- tmp2 - maxx

tmp2 <- exp(tmp2)

parm$clust$post.k <- tmp2
parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)

########################################################################
# store current state
old.parm <- parm

########################################################################

new.c.k <- sample(0:(parm$clust$G+1), size=1, replace=TRUE, prob=parm$clust$post.k)
parm$clust$c.v[k] <- new.c.k
new.flag <- new.c.k == (parm$clust$G+1)

#######################

count.0 <- sum(new.c.k==0)
parm$clust$C.m0 <- parm$clust$C.m0 + count.0

if ((count.0 == 0)&(!new.flag))
{parm$clust$C.m.vec[new.c.k] <- parm$clust$C.m.vec[new.c.k] + 1
}

if (new.flag)
{
  ###generate the latent vector first, condition on the single kth column
  cand.s.v.k <- array(,length(x.mt))
  
  for (tt in 1:length(x.mt))
  {
    tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)
    tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
    
    
    tmp.prob.v <- tmp.lik.v*P.aux
    prob.gen.v <- tmp.prob.v/sum(tmp.prob.v)
    cand.s.v.k[tt]<-sample(0:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v) # HOT
  } # end for loop
  
  parm$cand$s.v.k <- cand.s.v.k
  
  parm$cand$n.vec.k <- array(,parm$clust$K)
  for (gg in 1:parm$clust$K)
  {parm$cand$n.vec.k[gg] <- sum(cand.s.v.k==gg)
  }
  parm$cand$n0.k <- sum(cand.s.v.k==0)
  
  ##################
  parm$clust$G <- parm$clust$G + 1
  
  parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)
  
  parm$clust$beta.v <- c(parm$clust$beta.v, 0)
  parm$clust$gamma.v <- c(parm$clust$gamma.v,0)
  
  parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)
  
  parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)
  
  parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
  parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k
  
  parm$N <- sum(parm$clust$n.vec) + parm$clust$n0
  
  tmp.a.v <- array(,parm$n2)
  s.G.v <- parm$cand$s.v.k
  indxx <- s.G.v==0
  tmp.a.v[indxx] <- 0
  tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]]
  #
  parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
  # parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v) # HOT
  # 
  # if (parm$tBB_flag)
  # {
  # 
  #   parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
  # 
  # }
  
} # end  if (new.flag)


list(parm, new.flag)
}

###########################################################

row_PDP_fn.fast_col <- function(cc, parm, data)
{
  k <- parm$k <- parm$clust$col.nbhd.k[[cc]]
  I.k <- parm$clust$col.nbhd[[cc]]
  
  err <- row_PDP_fn.consistency.check(parm)
  if (err > 0)
  {stop(paste("FAST - 0: failed consistency check: err=",err))
  }
  
  # store so that we can revert to this state if MH propsal is rejected
  init.cc.parm <- parm
  
  old.c.k <- parm$clust$c.v[I.k]
  
  # Note to MS: I've ignored the branching conditions
  #             in "elementwise_DP.functions.R" based on computeMode$useR
  
  
  parm$clust$C.m.vec.k <- array(,parm$clust$G)
  for (gg in 1:parm$clust$G) {
    parm$clust$C.m.vec.k[gg] <- sum(old.c.k == gg)
    
    parm$clust$C.m0.k <- sum(old.c.k == 0)
  }
  
  
  parm$clust$C.m.vec.k.comp <- parm$clust$C.m.vec - parm$clust$C.m.vec.k
  parm$clust$C.m0.k.comp <- parm$clust$C.m0 - parm$clust$C.m0.k
  
  x.mt <- matrix(parm$X[,k], ncol=1)
  
  
  # intercept cluster or any existing cluster
  L.v <- sapply(0:parm$clust$G, row_PDP_fn.log.lik, x.mt, parm) # HOT
  
  emptied.indx <- which(parm$clust$C.m.vec.k.comp==0)
  new.G <- parm$clust$G - length(emptied.indx)
  
  if (length(emptied.indx) >0)
  {
    
    new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
    new.n.vec <- array(,parm$clust$K)
    for (pp in 1:parm$clust$K) {
      new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
    }
    new.n0 <- sum(new.s.mt == 0) # HOT 3
    
    emptied.s.indx <- which(new.n.vec == 0)
    new.K <- parm$clust$K - length(emptied.s.indx)
  }
  
  if (length(emptied.indx) ==0)
  {
    new.s.mt <- parm$clust$s.mt
    new.n.vec <- parm$clust$n.vec
    emptied.s.indx <- which(new.n.vec==0)
    new.K <- parm$clust$K
    new.n0 <- parm$clust$n0
  }
  
  ## generate auxilliary P vector
  
  tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  tmp.alpha <- tmp.M+new.n.vec
  tmp.alpha[emptied.s.indx] <- 0
  P.aux <- rgamma(parm$clust$K+1,c(0,tmp.alpha),1)
  # rgamma(parm$clust$K+1,c(parm$clust$M0+new.n0,tmp.alpha),1)
  P.aux <- P.aux/sum(P.aux)
  
  
  ## marginal likelihood of new cluster
  marg.log.lik.v <- cand.s.v.k <- array(,length(x.mt))
  for (tt in 1:length(x.mt))
  {
    tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) # HOT 3
    tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
    
    
    tmp.marg.v <- tmp.lik.v*P.aux
    marg.log.lik.v[tt] <- log(sum(tmp.marg.v))
    cand.s.v.k[tt]<-sample(0:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v)
  }
  
  parm$cand$s.v.k <- cand.s.v.k
  marg.log.lik <- sum(marg.log.lik.v)
  
  L.v <- c(L.v, marg.log.lik)
  
  ##
  log.prior.v <- array(NA, (2+parm$clust$G))
  
  spread.mass <- 0 #(parm$b1+new.G*parm$d)/(1+length(emptied.indx))
  
  if (length(emptied.indx) >0)
  {
    log.prior.v[-(emptied.indx+1)] <- log(c(0, (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
    # log(c((parm$b0+parm$clust$C.m0), (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
    log.prior.v[emptied.indx+1] <- log(spread.mass)
  }
  
  if (length(emptied.indx) ==0)
  {log.prior.v <- log(c(0, (parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
  # log(c((parm$b0+parm$clust$C.m0.k.comp), (parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
  }
  
  tmp2 <- log.prior.v + L.v
  maxx <- max(tmp2)
  tmp2 <- tmp2 - maxx
  
  tmp2 <- exp(tmp2)
  
  parm$clust$post.k <- tmp2
  parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)
  
  ########################################################################
  # store current state
  old.parm <- parm
  
  ########################################################################
  
  new.c.k <- sample(0:(parm$clust$G+1), size=length(I.k), replace=TRUE, prob=parm$clust$post.k)
  
  exit <- (sum(new.c.k != old.c.k)==0)
  flip <- TRUE
  new.flag <- FALSE
  
  if (!exit) # CONTINUE W/O EXITING FUNCTION, accepting M-H sampling
  {
    
    parm$clust$c.v[I.k] <- new.c.k
    
    new.count <- sum(new.c.k == (parm$clust$G+1))
    
    #######################
    
    new.prop <- 0
    
    count.0 <- sum(new.c.k==0)
    parm$clust$C.m0 <- parm$clust$C.m0.k.comp + count.0
    if (count.0 >0)
    {new.prop <- new.prop + log(parm$clust$post.k[1])*count.0
    }
    
    for (gg in 1:parm$clust$G)
    {count.gg <- sum(new.c.k==gg)
    parm$clust$C.m.vec[gg] <- parm$clust$C.m.vec.k.comp[gg] + count.gg
    
    if (count.gg > 0)
    {new.prop <- new.prop + log(parm$clust$post.k[gg+1])*count.gg
    }
    }
    
    
    if (new.count > 0)
    {
      parm$clust$G <- parm$clust$G + 1
      
      parm$clust$C.m.vec <- c(parm$clust$C.m.vec, new.count)
      
      new.prop <- new.prop + log(parm$clust$post.k[parm$clust$G+1])*new.count
      
      parm$cand$n.vec.k <- array(,parm$clust$K)
      for (ss in 1:parm$clust$K)
      {parm$cand$n.vec.k[ss] <- sum(parm$cand$s.v.k==ss)
      }
      parm$cand$n0.k <- sum(parm$cand.s.v.k==0)
      
      ##################
      
      
      parm$clust$beta.v <- c(parm$clust$beta.v, 0)
      parm$clust$gamma.v <- c(parm$clust$gamma.v,0)
      
      parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)
      
      parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$n2)
      
      parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
      parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k
      
      parm$N <- sum(parm$clust$n.vec) + parm$clust$n0
      
      tmp.a.v <- array(,parm$n2)
      s.G.v <- parm$cand$s.v.k
      indxx <- s.G.v==0
      tmp.a.v[indxx] <- 0
      tmp.a.v[!indxx] <- parm$clust$phi.v[s.G.v[!indxx]]
      #
      parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
      # parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)
      # 
      # if (parm$tBB_flag)
      # {
      #   parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
      # 
      # }
      
    } # end  if (new.count > 0)
    
    ######################
    
    # sum(parm$clust$C.m.vec) + parm$clust$C.m0 == parm$p
    
    
    
    ##########################################
    ##########################################
    ##### Computing proposal prob of reverse move
    ##########################################
    ##########################################
    
    old.prop <- 0
    
    for (gg in 0:old.parm$clust$G)
    {flag.gg <- (old.c.k==gg)
    count.gg <- sum(flag.gg)
    
    if (count.gg > 0)
    {old.prop <- old.prop + log(old.parm$clust$post.k[gg+1])*count.gg
    }
    }
    
    rho.prop <- new.prop - old.prop
    
    #######################################################
    #######################################################
    ########## computing true log-ratio: 2015
    #######################################################
    #######################################################
    
    # formula on page 1 of 07/27/15 notes
    
    # need to ensure there are no empty clusters in either
    # init.cc.parm or parm, otherwise likelihood formula
    # doesn't work (lgamma of negative values)
    
    tmp.new.parm <- parm
    indx.new <- parm$clust$C.m.vec > 0
    tmp.new.parm$clust$G <- sum(indx.new)
    tmp.new.parm$clust$C.m.vec <- parm$clust$C.m.vec[indx.new]
    #
    tmp.old.parm <- init.cc.parm
    indx.old <- init.cc.parm$clust$C.m.vec > 0
    tmp.old.parm$clust$G <- sum(indx.old)
    tmp.old.parm$clust$C.m.vec <- init.cc.parm$clust$C.m.vec[indx.old]
    
    rho.tru <- fn.d(d=parm$d, tmp.new.parm) - fn.d(d=parm$d, tmp.old.parm)
    
    new.log.lik <- 0
    
    for (gg in new.c.k)
    {indx.gg <- new.c.k==gg
    x_gg.mt <- matrix(parm$X[,I.k[indx.gg]], ncol=sum(indx.gg))
    new.log.lik <- new.log.lik + sum(row_PDP_fn.log.lik(gg, x.mt=x_gg.mt, parm))
    }
    
    old.log.lik <- 0
    for (gg in old.c.k)
    {indx.gg <- old.c.k==gg
    x_gg.mt <- matrix(parm$X[,I.k[indx.gg]], ncol=sum(indx.gg))
    old.log.lik <- old.log.lik + sum(row_PDP_fn.log.lik(gg, x.mt=x_gg.mt, old.parm))
    }
    
    rho.tru <- rho.tru + new.log.lik - old.log.lik
    
    ########## toss a coin #################
    
    prob <- exp(min((rho.tru-rho.prop),0))
    
    flip<-as.logical(rbinom(n=1,size=1,prob=prob))
    
    if (!flip) {parm <- init.cc.parm}
    
    
  } # end BIG if (!exit) loop
  
  
  list(parm, new.flag, exit, flip)
}


#####################################


row_PDP_fn.drop <- function(parm, computeMode)
{
  
  ##########################################
  ## Drop empty clusters:
  ## (i)  Move empty clusters to end by relabeling clusters
  ## (ii) Set parm$clust$G equal to number of non-empty clusters
  ## (ii) Retain only clusters  1,...,parm$clust$G
  #########################################
  
  parm$clust$G <- sum(parm$clust$C.m.vec>0)
  num.dropped <- sum(parm$clust$C.m.vec==0)
  
  if (parm$clust$G > 0)
  {
    if (num.dropped > 0)
    {
      for (rr in 1:num.dropped)
      {
        old.label <-  min(which(parm$clust$C.m.vec==0))
        new.label <- max(which(parm$clust$C.m.vec>0))
        stopp <-  max(which(parm$clust$C.m.vec>0)) == parm$clust$G
        if (stopp)
        {break
        }
        parm <- row_PDP_fn.swap.clusters(parm, g1 = new.label, g2 = old.label)
      }
    }
    
    ##########
    
    keep <- 1:parm$clust$G
    parm <- row_PDP_fn.clip.clusters(parm, keep)
    
    ###########
    
    # parm$clust$K does not change (possibly some empty elementwise clusters)
    
    parm$N <- parm$n2 * parm$clust$G
    parm$clust$s.v <- as.vector(parm$clust$s.mt)
    
    parm$clust$n0 <- sum(parm$clust$s.v==0)
    parm$clust$n.vec <- array(,parm$clust$K)
    for (ss in 1:parm$clust$K)
    {parm$clust$n.vec[ss] <- sum(parm$clust$s.v==ss)  # HOT
    }
    
  }
  parm
}



row_PDP_fn.main <- function(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size)
{
  
  parm$n2 <- parm$clust$GG
  p <- parm$p
  
  
  ##########################
  # compute delta-neighborhoods
  #########################
  
  # SG: prob.compute.col.nbhd is the probability of computing
  # the neighborhoods, and col.frac.probes is the fraction of neighborhoods updated
  
  col_flip <- as.logical(rbinom(n=1,size=1,prob=prob.compute.col.nbhd))
  if (is.null(parm$clust$col.nbhd.k) | col_flip){
    parm <- row_PDP_fn.post.prob.and.delta(parm, max.col.nbhd.size, col.frac.probes)
  }
  
  
  if (col.frac.probes < 1)
  {num_nbhds <- max(1,round(col.frac.probes*length(parm$clust$col.nbhd.k)))
  parm$subset_nbhd.indx <- sort(sample(1:length(parm$clust$col.nbhd.k), size=num_nbhds))
  }
  
  if (col.frac.probes == 1)
  {parm$subset_nbhd.indx <- 1:length(parm$clust$col.nbhd.k)
  }
  
  new.flag.v <- NULL
  col.mh.flip.v <- col.mh.exit.v <- NULL
  
  for (cc in parm$subset_nbhd.indx)
  {previous.parm <- parm
  
  parm$k <- parm$clust$col.nbhd.k[[cc]]
  
  if (length(parm$clust$col.nbhd[[cc]])==1)
  { tmp <- row_PDP_fn.gibbs(k=parm$k, parm, data)
  parm <- tmp[[1]]
  new.flag.v <- c(new.flag.v, tmp[[2]])
  }
  
  if (length(parm$clust$col.nbhd[[cc]])>1)
  {
    tmp <- row_PDP_fn.fast_col(cc, parm, data)
    parm <- tmp[[1]]
    new.flag.v <- c(new.flag.v, tmp[[2]])
    col.mh.exit.v <- c(col.mh.exit.v, tmp[[3]])
    col.mh.flip.v <- c(col.mh.flip.v, tmp[[4]])
  }
  
  } # end for loop
  
  
  err <- row_PDP_fn.consistency.check(parm)
  if (err > 0)
  {stop(paste("LOOP: failed consistency check: err=",err))
  }
  
  parm$clust$col.new.flag <- mean(new.flag.v)
  if (!is.null(col.mh.flip.v))
  {parm$clust$col.mh.flip <- mean(col.mh.flip.v)
  parm$clust$col.mh.exit <- mean(col.mh.exit.v)
  }
  else
  {parm$clust$col.mh.flip <- 1
  parm$clust$col.mh.exit <- 1
  }
  
  ##########################################
  ## Drop empty group clusters:
  ## (i)  Move empty clusters to end by relabeling clusters
  ## (ii) Set parm$clust$G equal to number of non-empty clusters
  ## (ii) Retain only clusters  1,...,parm$clust$G
  #########################################
  
  parm <- row_PDP_fn.drop(parm)
  
  # now drop empty elementwise clusters
  
  parm <- row_element_fn.drop(parm)
  
  err <- row_PDP_fn.compact.consistency.check(parm)
  if (err > 0)
  {stop(paste("MAIN FUNCTION END: failed consistency check: err=",err))
  }
  
  parm
}



