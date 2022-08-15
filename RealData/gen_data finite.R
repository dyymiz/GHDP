gen.X.finite <- function(n, prop.X.miss, true_parm) {
  
  data <- NULL
  data$X <-  rep(list(array(,c(n,p))),t)
  data$Xorder <- rep(list(array(,c(n,p))),t)
  true_parm$theta.mt<- array(,c(n,p,t))
  # pr <- runif(m,0,10)
  pr <- c(rep(1,(m-1)),6*m)
  true_parm$pr <- pr/sum(pr)

  true_parm$clust$r.v <- sample(1:m,prob=pr,size=n,replace=TRUE)
  r <- true_parm$clust$R.m.vec <- table(true_parm$clust$r.v)
  
  # for (tt in 1:t) {
  #   for (mm in 1:m) {
  #     for (xx in 1:p){
  #     data$X[[tt]][,xx] <- rnorm(n=n,mean=true_parm$theta[components,xx,tt],sd=true_parm$tau)
  #     }
  #   }
  # }
  
  for(tt in 1:t){
    true_parm$theta.mt[,,tt] <- apply(true_parm$theta[,,tt],2, function(x) rep(x,times=r))
   for (xx in 1:p)
    {
     data$X[[tt]][,xx] <- rnorm(n=n,mean=true_parm$theta.mt[,xx,tt],sd=true_parm$tau)
    }
  }
  
  tmp.mat <- array(0,c(n,n))
  
  for (jj in unique(true_parm$clust$r.v))
  {indx.jj <- which(true_parm$clust$r.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  true_parm$clust$sub.nbhd.matrix <- tmp.mat
  
  ###########################################
  # missing X values
  ###########################################
  
  data$num.X.miss <- round(prop.X.miss*n*p)
  
  if (data$num.X.miss>0)
  {data$X.missing.x <- sample(1:n,size=data$num.X.miss, replace=TRUE)
  data$X.missing.y <- sample(1:p,size=data$num.X.miss, replace=TRUE)
  }
  
  # works even if some (data$X.missing.x,data$X.missing.y) are tied by chance
  for(tt in 1:t){
   for (cc in 1:data$num.X.miss)
    {
     data$X[[tt]][data$X.missing.x[cc],data$X.missing.y[cc]] <- NA
    }
  }
  
  data$K.max <- round(true_parm$clust$alpha_0*log(p))
  data$G.max <- round(p/2)
  # dimension os data$X[[tt]] doesn't match with n2 
  
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

