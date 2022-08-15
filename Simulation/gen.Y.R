
fn.comb.s <- function(All.Stuff){
  
  All.Stuff.s <- NULL
  All.Stuff.s$parm <- All.Stuff$parm.m

  tmp.A <- do.call("rbind", All.Stuff$est_c.v.c)
  tmp_c.v <- c.tmp <- tmp.A
  c.tmp[,1] <- lapply(1:s, function(x) do.call(max,tmp.A[x,1]))
  for (i in 2:t) {
    tmp_c.v[,i] <- lapply((i-1):s,function(x) do.call("+", c.tmp[x,]))
    c.tmp[,i] <- lapply(1:s, function(x) do.call(max,tmp.A[x,1]))
  }
  tmp.A <- tmp_c.v
  
  All.Stuff.s$parm$est_c.v <- lapply(1:s, function(x) do.call("c",tmp.A[x,]))
  All.Stuff.s$parm$clust$G <- lapply(1:s, function(x) max(All.Stuff.s$parm$est_c.v[[x]]))

  All.Stuff.s
}

fn.slice.s <- function(All.Stuff.s, mm){
  
  All.Stuff.m <- NULL
  All.Stuff.m$parm$clust$G <- All.Stuff.s$parm$clust$G[[mm]]
  # All.Stuff.m$parm$clust$A.mt <- All.Stuff.s$parm$clust$A.mt[[mm]]
  # All.Stuff.m$est_c.v <- All.Stuff.s$est_c.v[[mm]]
  All.Stuff.m$parm$est_c.v <- All.Stuff.s$parm$est_c.v[[mm]]
  All.Stuff.m$parm$X <- All.Stuff.s$parm[[mm]]$X
  All.Stuff.m
  
}

fn_gen.Y.summary <- function(parm)
{	
  parm$clust$c.v <- parm$est_c.v
  parm$clust$G <- max(parm$est_c.v)
  
  parm$clust$C.m.vec <- array(,parm$clust$G)
  
  for (g in 1:parm$clust$G)
  {I.g <- (parm$clust$c.v==g)
  parm$clust$C.m.vec[g] <- m.g <- sum(I.g)
  }
  
  parm$clust$C.m0 <- p - sum(parm$clust$C.m.vec)
  
  #################################
  
  parm
}


fn.gen.Y <- function(All.Stuff, data, beta, cens.prop, num.pred)
{
  #########################################
  # generating the responses Y
  ########################################
  
  parm <- All.Stuff$parm
  parm <- fn_gen.Y.summary(parm)
  
  # select predictors' cluster that have small sizes (=> typically, larger correlations)
  true.pred.clusters.v <- sample(1:parm$clust$G, prob=1/parm$clust$C.m.vec^2, size=num.pred)
  
  # pick true predictors from member covariates of those predictor cluster
  true.pred.v <- array(,num.pred)
  
  for (xx in 1:num.pred)
  {indx.xx <- which(parm$clust$c.v==true.pred.clusters.v[xx])
  if (length(indx.xx)>1)
  {true.pred.v[xx] <- sample(indx.xx, size=1)
  }
  if (length(indx.xx)==1)
  {true.pred.v[xx] <- indx.xx
  }
  }
  
  true.pred.v <- sort(true.pred.v)
  
  trueY = NULL
  trueY$true.pred.v <- true.pred.v 
  trueY$num.pred <- num.pred
  trueY$true.pred.clusters.v <- sort(true.pred.clusters.v)
  
  trueY$beta.v <- rep(0,p)
  trueY$beta.v[true.pred.v] <- beta
  trueY$beta0 <- 0
  
  # note: parm$X is just data$X with any missing values imputed
  # trueY$mean.v <- trueY$beta0 + as.vector(parm$X[,true.pred.v] %*% trueY$beta.v[true.pred.v])
  # trueY$pa <- parm
  
  rm(parm)
  
  # change this to real data, and non.missing.indx, missing.indx, and delta
  # trueY$Y <- log(rgamma(n=n2, shape=1, scale=exp(trueY$mean.v)))
  trueY$Y <- log(clin[,"days_to_death"])[n_sub]
  # all uncensored delta
  # trueY$delta <- rep(1, n2)
  trueY$delta <- clin[,"censor"][n_sub]
  
  # if (cens.prop > 0)
  # 	{
  # 	# first round(n2*cens.prop) subjects censored
  # 	trueY$delta[1:round(n2*cens.prop)] <- 0 # indicator of censoring
  # 
  # 	for (i in 1:round(n2*cens.prop))
  # 		{avg <- trueY$mean.v[i]
  # 		prob <- pgamma(exp(trueY$Y[i]), shape=1, scale=exp(avg))
  # 
  # 		flag1 <- (prob > 1e-5)
  # 		flag2 <- (prob < (1- 1e-5))
  # 		 
  # 		if ((flag1)&(flag2))
  # 			{u <- runif(n=1, min=0, max=prob)
  # 			tmp.y <- log(qgamma(u, shape=1, scale=exp(avg)))
  # 			}
  # 		if ((!flag1)|(!flag2))
  # 			{stop("Error in response generation")
  # 			}
  # 		trueY$Y[i] <- tmp.y
  # 		}
  # 	}
  
  ###########################################
  # random split of 100 X test.prop % missing
  ###########################################
  
  indx <- sort(sample(1:n2, size=n1))
  data$missing.indx <- indx
  data$non.missing.indx <- setdiff(1:n2, data$missing.indx)
  
  data$Y <- trueY$Y[data$non.missing.indx]
  data$delta <- trueY$delta[data$non.missing.indx]
  data$true <- NULL
  data$true$Y <- data$trueY$Y <- trueY$Y
  data$true$delta <- data$trueY$delta <- trueY$delta
  
  if (n1==0) 
  {data$Y <- trueY$Y
  data$delta <- trueY$delta
  } 
  
  data$Y <- trueY$Y[data$non.missing.indx]
  data$delta <- trueY$delta[data$non.missing.indx]
  data$true <- NULL
  data$true$Y <- data$trueY$Y <- trueY$Y
  data$true$delta <- data$trueY$delta <- trueY$delta
  
  list(trueY, data)
  
}


