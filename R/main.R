rm(list=ls())

library(MASS)
library(mvtnorm)
library(MCMCpack)
library(expm)
library(stringr)
library(abind)
library(fields)
library(plyr)
library(survival)
library(mboost)
library(coda)


options(error=recover)
options(warn=2)

source("sourceAll.R")

path = getwd()
sourceDir(path, except=c("main.R","main_TCGA.R","concordance.r","interpret.R", "subtype_summary.R"))
load("LungCancerDataLUSCLUADClinical.RData")


t=2
s=1
m=3
w=500

n <- dim(data[[1]])[1]
p <- array(c(dim(data[[1]])[2], dim(data[[2]])[2]),c(s,t))
r <- n/s

prop.X.miss=0
tau = 0.2
mu.2=0
tau.2=1.25

do.response = TRUE
beta.v <- 1
cens.prop.v <- .2
y.reps = 1 
test.prop <- 0
num.pred <- 20
num.knots <- 1
do.gamma.prob <- 1/20

true_parm <- data1 <- NULL
true_parm$clust$alpha_0 <- true_parm$clust$alpha_1 <- 22
n.burn = 500
n.reps = 1000
max.row.nbhd.size = round(.1*25*250^.5) 
max.col.nbhd.size = round(.05*250) 
row.frac.probes = 0.05
col.frac.probes = .1
prob.compute.col.nbhd=.2
tBB_flag=FALSE

data <- lapply(1:t, function(x) t(apply(data[[x]], 1, scale)))
lapply(1:t, function(x) image.plot(t(data[[x]])))

data1$X[[1]] <- data[[1]]
data1$X[[2]] <- data[[2]]
data <- data1

lapply(1:t, function(x) image.plot(t(data1$X[[x]])))

time.elapse <- Sys.time()

All.Stuff <- fn.mcmc(true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, 
              row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, tBB_flag)

p <- sum(unlist(lapply(1:t, function(x) dim(data1$X[[x]])[2])))

All.Stuff.s <- fn.comb.s(All.Stuff)
time.elapse <- Sys.time()

if (do.response){
  
  num.pred.mt <- train.r.sq.mt <-  train.mse.mt <-  train.abs.mt <-  mh.accept.mt  <- train.concord.mt <- test1.concord.mt <- test2.concord.mt <- test3.concord.mt <- test4.concord.mt <- test5.concord.mt <- test6.concord.mt <- test7.concord.mt <- prob.second.knot.mt  <- array(,c(s, y.reps,length(cens.prop.v),length(beta.v)))
  
  prob.knot.13.ar <-  array(,c(s, y.reps,length(cens.prop.v),length(beta.v),c(2,2)))
  
  fdr.mt <- fnr.mt <-  array(,c(s, y.reps,length(cens.prop.v),length(beta.v)))
  
  max.Z.mt <- change.gamma.mt <- array(,c(s, y.reps,length(cens.prop.v),length(beta.v)))
  
  cor1.missing.mt <- cor2.missing.mt <- cor.non_missing.mt <- r1.sq.missing.mt <- r2.sq.missing.mt <- r.sq.non_missing.mt <-  array(,c(s, y.reps,length(cens.prop.v),length(beta.v)))
  
  mse1.missing.mt <- mse2.missing.mt <- mse.non_missing.mt <- array(,c(s, y.reps,length(cens.prop.v),length(beta.v)))
  
  
  for(mm in 1:s){
    All.Stuff.m <- fn.slice.s(All.Stuff.s,mm)
    print(" ")
    print("****************************")
    print(" ")
    print(paste("s = ",mm,date(),"***********"))
    print(" ")
    print("****************************")
    print(" ")

    n2 = r[[mm]]
    n1 <- round(test.prop*n2)
    n0 = n2-n1
    
    n.burn <- 10000
    n.reps <- 8000
    
    ###################
    # define storage objects
    ###################
    
    train.mean.ar <- test.mean.ar <- array(,c(n2,y.reps, length(cens.prop.v),length(beta.v)))
    
    ###################
    # iterations
    ###################
    
    for (xxx in 1:y.reps)
      for (yyy in 1:length(cens.prop.v))
        for (zzz in 1:length(beta.v))
        {
          set.seed(1000 + 200 + 9*xxx)
          
          beta <- beta.v[zzz]
          cens.prop <- cens.prop.v[yyy]
          
          simY = fn.gen.Y(All.Stuff.m, data, beta, cens.prop, num.pred)
          
          trueY = simY[[1]]
          data = simY[[2]]
          
          rm(simY)
          
          
          ###################
          # given covariates and clusters, generate responses
          # and detect predictors
          ###################
          
          All.Stuff.Y <- fn_varsel.mcmc(All.Stuff.m, true, data, n.burn, n.reps, num.knots, do.gamma.prob)
          
          ########################
          
          prob.knot.13.ar[mm, xxx,yyy,zzz,,] <- All.Stuff.Y$gamma_13.mt
          
          ########################
          
          max.Z.mt[mm,xxx,yyy,zzz] <- mean(All.Stuff.Y$maxed.Z.flag.v)
          
          change.gamma.mt[mm,xxx,yyy,zzz] <- mean(All.Stuff.Y$change.gamma.v, na.rm=TRUE)
          
          indx.1 <- which(data$true$delta==1); indx.nm <- data$non.missing.indx
          indx10 <- intersect(which(data$true$delta==1), data$missing.indx)
          indx11 <- intersect(which(data$true$delta==1), data$non.missing.indx)
          
          cor1.missing.mt[mm,xxx,yyy,zzz] <- cor(data$true$Y[indx10], All.Stuff.Y$test.mean.v[indx10])
          cor2.missing.mt[mm,xxx,yyy,zzz] <- cor(data$true$Y[indx10], All.Stuff.Y$mean.v[indx10])
          cor.non_missing.mt[mm,xxx,yyy,zzz] <- cor(data$true$Y[indx11], All.Stuff.Y$mean.v[indx11])
          
          r1.sq.missing.mt[mm,xxx,yyy,zzz] <- 1 - mean((data$true$Y[indx10] - All.Stuff.Y$test.mean.v[indx10])^2) / var(data$true$Y[indx10])
          r2.sq.missing.mt[mm,xxx,yyy,zzz] <- 1 - mean((data$true$Y[indx10] - All.Stuff.Y$mean.v[indx10])^2) / var(data$true$Y[indx10])
          r.sq.non_missing.mt[mm,xxx,yyy,zzz] <-  1 - mean((data$true$Y[indx11] - All.Stuff.Y$mean.v[indx11])^2) / var(data$true$Y[indx11])
          
          mse1.missing.mt[mm,xxx,yyy,zzz] <- mean((data$true$Y[indx10] - All.Stuff.Y$test.mean.v[indx10])^2)
          mse2.missing.mt[mm,xxx,yyy,zzz] <- mean((data$true$Y[indx10] - All.Stuff.Y$mean.v[indx10])^2)
          mse.non_missing.mt[mm,xxx,yyy,zzz] <- mean((data$true$Y[indx11] - All.Stuff.Y$mean.v[indx11])^2)
          
          ########################
          
          num.pred.mt[mm, xxx,yyy,zzz] <- mean(All.Stuff.Y$num.pred.v)
          
          indx11 <- intersect(which(data$true$delta==1), data$non.missing.indx)
          # gain of using variable selection 
          # relative to no-predictor model FOR TRAINING CASES
          
          train.r.sq.mt[mm,xxx,yyy,zzz] <- 1 - mean((data$true$Y[indx11] - All.Stuff.Y$mean.v[indx11])^2) / var(data$true$Y[indx11]) 
          train.mse.mt[mm,xxx,yyy,zzz] <- mean((data$true$Y[indx11] - All.Stuff.Y$mean.v[indx11])^2)
          train.abs.mt[mm,xxx,yyy,zzz] <- mean(abs(data$true$Y[indx11] - All.Stuff.Y$mean.v[indx11]))
          
          ########################
          
          train.mean.ar[,xxx,yyy,zzz] <- All.Stuff.Y$mean.v
          test.mean.ar[,xxx,yyy,zzz] <- All.Stuff.Y$test.mean.v
          
          indx.set <- data$non.missing.indx	
          indx2.set <- data$non.missing.indx
          mean.v <- All.Stuff.Y$mean.v
          source("concordance.r")
          train.concord.mt[mm,xxx,yyy,zzz] <- error.rate
        
          
          print(" ")
          print("****************************")
          print("****************************")
          print(" ")
          print(paste("FINISHED"))
          print(" ")
          print("****************************")
          print("****************************")
          print(" ")
          
        }
  
  }
  
}# end massive if (do.response) loop

Sys.time()-time.elapse

B_fdr<-function(bf_vec, alpha=0.35){
  bfv = sort(bf_vec, decreasing = TRUE)
  cm = cumsum(1- bfv)
  tmp_alpha = max(which(cm < alpha))
  b_fdr = bfv[tmp_alpha]
  b_fdr
}
b_fdr = B_fdr(All.Stuff.Y$is.pred.v)
b_fdr

# selected probs
selected.clst <- which(All.Stuff.Y$is.pred.v >= 0.01)
sum(unlist(All.Stuff.s$parm$est_c.v) %in% selected.clst)
sort(All.Stuff.Y$is.pred.v[selected.clst])


# selected probes in each platform
clst.mRNA <- selected.clst[selected.clst <= max(All.Stuff$est_c.v.c[[1]][[1]])]
clst.meth <- selected.clst[selected.clst > max(All.Stuff$est_c.v.c[[1]][[1]])] - max(All.Stuff$est_c.v.c[[1]][[1]])
length(clst.mRNA)
length(clst.meth)

sum(All.Stuff$est_c.v.c[[1]][[1]] %in% clst.mRNA)
sum(All.Stuff$est_c.v.c[[1]][[2]] %in% clst.meth)

# select one prob randomly from every cluster
prob.mRNA <- prob.meth <- c()

for(i in 1:length(clst.mRNA)){
  mRNA_tmp <- which(All.Stuff$est_c.v.c[[1]][[1]]==clst.mRNA[i])
  
  if(length(mRNA_tmp) > 1){
    prob.mRNA[i] <- sample(x=mRNA_tmp, size=1)
  }else{prob.mRNA[i] = mRNA_tmp}
}

for(i in 1:length(clst.meth)){
  meth_tmp <- which(All.Stuff$est_c.v.c[[1]][[2]]==clst.meth[i])
  
  if(length(meth_tmp) > 1){
    prob.meth[i] <- sample(meth_tmp, 1)
  }else{prob.meth[i] = meth_tmp}
}

load("RealData.RData")
load("LungCancerDataLUSCLUADClinical.restart.2020.RData")

sub.genes <- subset(merged.data, Platform == "geneExp")
sub.methylation <- subset(merged.data, Platform == "methylation")

genes.RNA <- sub.genes[colnames(data[[1]]), 1]
genes.meth <- sub.methylation[colnames(data[[2]]), 1]

slt.RNA <- genes.RNA[prob.mRNA]
slt.meth <- genes.meth[prob.meth]
slt.RNA
slt.meth
intersect(slt.RNA,slt.meth)

all.genes <- list()
for(i in 1:length(clst.mRNA)){
  all.genes[[i]] <- c(genes.RNA[which(All.Stuff$est_c.v.c[[1]][[1]]==clst.mRNA[i])])
}

all.meth <- list()
for(i in 1:length(clst.meth)){
  all.meth[[i]] <- c(genes.meth[which(All.Stuff$est_c.v.c[[1]][[2]]==clst.meth[i])])
}

