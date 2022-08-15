SBP <- function(num_weights, alpha,true_parm) {
  
  betas <- rbeta(num_weights, 1, alpha)
  remaining_stick_lengths <- c(1, cumprod(1 - betas))[1:(num_weights-1)]
  true_parm$weights <- remaining_stick_lengths * betas[1:(num_weights-1)] 
  # truncated weights add up to be 1
  true_parm$weights <- c(true_parm$weights,1-sum(true_parm$weights))
  true_parm$phi <- rnorm(num_weights,mu.2, tau.2)
  # labels of unique phi_{k^*}
  true_parm$clust$t.v <- c(1:num_weights)
  
  true_parm
}

gen.clust <- function(s, w, p) {
  # m is # of row clusters
  # w is # of phi_{mst} mapped from phi_{s}
  true_parm <- NULL
  true_parm$clust <- NULL  # column cluster membership counts. Column cluster allocation variable

  # DP
  true_parm$clust$alpha_0 <- 5
  # PDP
  true_parm$clust$alpha_1 <- 5
  true_parm$s <- s
  # true_parm$d <- matrix(runif(s*t,0,1),nrow = s)
  true_parm$d <- matrix(c(0.2,0.25,0.30,0.35),nrow = s)

  true_parm$b0 <- 0
  true_parm$b1 <- 1
  
  true_parm$clust$n.vec <- 1
  true_parm$clust$n0 <- 0
  true_parm$clust$K <- 1
  true_parm$clust$M0 <- 0
  true_parm$clust$M <- true_parm$clust$alpha_0
  
  true_parm$clust$d.v <- array(, s*w)
  true_parm$clust$C.m.vec <- true_parm$phi.mt <- true_parm$clust$A.mt <- true_parm$clust$s.mt <- true_parm$clust$s.v <- true_parm$clust$c.v <- true_parm$clust$G <-rep(list(list()),t)
  true_parm$clust$col.nbhd.matrix <- rep(list(list()),t)
  true_parm.s <- true_parm$clust$r.v <- true_parm$clust$R.n.vec <- true_parm$clust$row.nbhd.matrix <- list()
  true_parm <- SBP(2*w,true_parm$clust$alpha_0,true_parm)

  # s cancer types
  for(ss in 1:s){
    for(tt in 1:t){
      # betas <- rbeta(1, 1-true_parm$d[ss,tt], true_parm$clust$alpha_1+1*true_parm$d[ss,tt])
      betas <- rbeta(1, 1, true_parm$clust$alpha_1)
      
      for(i in 2:w){
       # betas <- c(betas, rbeta(1, 1-true_parm$d[ss,tt], true_parm$clust$alpha_1+i*true_parm$d[ss,tt]))
        betas <- c(betas, rbeta(1, 1, true_parm$clust$alpha_1))
        
        }
      true_parm$clust$d.v <- matrix(sample(true_parm$clust$t.v,s*w,replace=TRUE, prob=true_parm$weights),nrow=s,byrow = TRUE)
      
      remaining_stick_lengths <- c(1, cumprod(1 - betas))[1:(w-1)]
      weights <- remaining_stick_lengths * betas[1:(w-1)] 
      weights <- c(weights,1-sum(weights))
      # map phi_{s} to phi_{mst}
      true_parm$phi.mt[[ss]][[tt]] <- true_parm$phi[true_parm$clust$d.v[ss,]]

      c.v <- array(,p[ss,tt])
      c.v[1] <- 1
      C.m.vec <- 1
      G <- 1
      
      for (xx in 2:p[ss,tt])
      {prob.v <- C.m.vec - true_parm$d[ss,tt]
      prob.v <- c(prob.v, (true_parm$b1 + G*true_parm$d[ss,tt]))
      new.c <- sample(1:(G+1), size=1, prob=prob.v)
      
      c.v[xx] <- new.c
      new.flag <- (new.c > G)
      
      if (new.flag)
      {G <- G + 1
      C.m.vec <- c(C.m.vec, 1)
      }
      if (!new.flag)
      {C.m.vec[new.c] <- C.m.vec[new.c] + 1
      }
      }
      
      s.v <- sample(true_parm$clust$d.v[ss,],m*G,replace = TRUE, prob = weights)
      s.mt <- matrix(s.v, nrow = m)
      
      true_parm$clust$c.v[[ss]][[tt]] <- c.v
      true_parm$clust$G[[ss]][[tt]] <- G
      true_parm$clust$C.m.vec[[ss]][[tt]] <- C.m.vec
      true_parm$clust$s.v[[ss]][[tt]] <- s.v
      true_parm$clust$s.mt[[ss]][[tt]] <- s.mt
      true_parm$clust$A.mt[[ss]][[tt]] <- matrix(true_parm$phi[s.mt],nrow = m)
      
      tmp.mat <- array(0,c(p[ss,tt],p[ss,tt]))

      for (jj in unique(true_parm$clust$c.v[[ss]][[tt]]))
      {indx.jj <- which(true_parm$clust$c.v[[ss]][[tt]]==jj)
      tmp.mat[indx.jj,indx.jj] <- 1
      }

      true_parm$clust$col.nbhd.matrix[[ss]][[tt]] <- tmp.mat
      
    } # end for t
    

    
    pr <- rep(1,m)
    # pr <- c(rep(1,(m-1)),m^2)
    true_parm$pr <- pr/sum(pr)
    
    true_parm$clust$r.v[[ss]] <- sample(1:m,prob=pr,size=r[ss],replace=TRUE)
    true_parm$clust$R.m.vec[[ss]] <- table(true_parm$clust$r.v[[ss]])
    
    tmp.mat <- array(0,c(r[ss],r[ss]))

    for (jj in unique(true_parm$clust$r.v[[ss]]))
    {indx.jj <- which(true_parm$clust$r.v[[ss]]==jj)
    tmp.mat[indx.jj,indx.jj] <- 1
    }

    true_parm$clust$row.nbhd.matrix[[ss]] <- tmp.mat
    
    # true_parm.s[[ss]] <- true_parm
  
    } # end for s
  
  ## neighborhood taxicab distances for column clusters
  
  
  ###########
  
  true_parm
}

