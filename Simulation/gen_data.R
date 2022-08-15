gen.X <- function(n, prop.X.miss, true_parm) {
  
  data <- NULL
  data$X <- data$true_theta <- rep(list(list()),t)

  for (ss in 1:s) {
    # true_parm <- true_parm.s[[ss]]
    true_parm$tau <- tau
    tmp <- list()
    R <- true_parm$clust$R.m.vec[[ss]]
    tmp[[1]] <- 1:cumsum(R)[1]
    for(rr in 1:(m-1)){
      tmp[[rr+1]] <- (cumsum(R)[rr]+1):cumsum(R)[rr+1]
    }
  for(tt in 1:t){
    
    X <- array(,c(r[ss],p[ss,tt]))
    theta.mt <- true_parm$clust$A.mt[[ss]][[tt]][true_parm$clust$r.v[[ss]],true_parm$clust$c.v[[ss]][[tt]]]
    
    for (j in 1:p[ss,tt]) {
      X[,j] <- rnorm(n=r[ss],mean=theta.mt[,j], sd=true_parm$tau)
      data$num.X.miss <- round(prop.X.miss*r[ss]*p[ss,tt])
    }
     if (data$num.X.miss>0)
     {data$X.missing.x <- sample(1:r[ss],size=data$num.X.miss, replace=TRUE)
     data$X.missing.y <- sample(1:p[ss,tt],size=data$num.X.miss, replace=TRUE)
     }
       for (cc in 1:data$num.X.miss)
       {
         X[data$X.missing.x[cc],data$X.missing.y[cc]] <- NA
       }
    data$X[[ss]][[tt]] <- X
    data$true_theta[[ss]][[tt]] <- theta.mt
    } # for t
    

  }# end for s
  
  if(s==1){
    data$X <- unlist(data$X, recursive = FALSE)
    data$true_theta <- unlist(data$true_theta, recursive = FALSE)
  }else{
  tmp.x=do.call("rbind", data$X)
  data$X <- lapply(1:s, function(x) do.call("rbind",tmp.x[,x]))
  
  tmp.x=do.call("rbind", data$true_theta)
  data$true_theta <- lapply(1:s, function(x) do.call("rbind",tmp.x[,x]))
  }
  
  p <- mean(p[1,])
  data$K.max <- round(true_parm$clust$alpha_0*log(p))
  data$G.max <- round(p/2)
  ############
  
  true <- NULL
  true$a.R <- true_parm$clust$M
  true$b0 <- 2.2
  true$b1 <- true_parm$b1
  
  #########################################
  # generating the R- and C- clusters
  ########################################
  
  true$shift <- 1e-4
  
  return(list(data = data, true = true, true_parm=true_parm))
}


fn.data.slice <- function(data,parm,true_parm,mm,tt,prop.X.miss=0){

  tmp_parm <- parm

  tmp_parm$clust$phi.v <- parm$clust$phi.list[[tt]]

    # relabeling c.v, could be a problem
    # tmp_parm$clust$c.v <- as.numeric(factor(rank(parm$clust$c.c.v[mm,,tt])))
    if(sum(parm$clust$c.v[[tt]]==0)>0){
      tmp_parm$clust$c.v <- as.numeric(factor(rank(parm$clust$c.v[[tt]])))-1
    }else{tmp_parm$clust$c.v <- as.numeric(factor(rank(parm$clust$c.v[[tt]])))}

  
  tmp_parm$clust$s.s.v <- parm$clust$s.v[[tt]]
  if(sum(tmp_parm$clust$s.s.v==0)>0){
    tmp_parm$clust$s.v <- as.numeric(factor(rank(tmp_parm$clust$s.s.v)))-1
  }else{tmp_parm$clust$s.v <- as.numeric(factor(rank(tmp_parm$clust$s.s.v)))}
  
  tmp_parm$clust$s.mt <- matrix(tmp_parm$clust$s.v, ncol = parm$clust$G[[tt]])
  tmp_parm$clust$s.s.mt <- parm$clust$s.mt[[tt]]
  tmp_parm$clust$G <- parm$clust$G[[tt]]
  # tmp_parm$clust$K <- max(tmp_parm$clust$s.v[[tt]])
  tmp_parm$clust$K <- parm$clust$K[[tt]]
  tmp_parm$clust$A.mt <- parm$clust$A.mt[[tt]]
  # check it 
  # tmp_parm$clust$n.vec <- parm$clust$n.vec[[tt]]
  tmp_parm$clust$n.vec <- tabulate(tmp_parm$clust$s.v)
  tmp_parm$clust$n.vec <- tmp_parm$clust$n.vec[tmp_parm$clust$n.vec!=0]
  tmp_parm$clust$C.m.vec <- array(,tmp_parm$clust$G)
  
  # labeling problem, may need to reorder
  for (g in 1:tmp_parm$clust$G)
  {I.g <- (tmp_parm$clust$c.v==g)
  tmp_parm$clust$C.m.vec[g] <- sum(I.g)
  }
  
      data1 <- NULL
      data1$X <- data
      n <- dim(data)[1]
      tmp_parm$p <- p <- dim(data)[2]
      tmp_parm$n2 <- n
      n1 <- 0
      
      tmp_parm$clust$C.m0 <- tmp_parm$p - sum(tmp_parm$clust$C.m.vec)
      
      data1$num.X.miss <- round(prop.X.miss*n*p)
      data1$missing.indx <- NULL
      data1$non.missing.indx <- 1:n
      
      data1$K.max <- round(true_parm$clust$alpha_0*log(tmp_parm$p))
      data1$G.max <- round(p/2)
      
      ###########################################
      # dummy responses
      ###########################################
      
      data1$Y <- rep(0,tmp_parm$n2)
      data1$delta <- rep(0,tmp_parm$n2)
      data1$true <- NULL
      data1$true$Y <- data1$Y
      data1$true$delta <- data1$delta

      # tmp_parm = tmp_parm, rewrite everything, naming problem
  return(list(data1=data1, tmp_parm=tmp_parm))
  
}
