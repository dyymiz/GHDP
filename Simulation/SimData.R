
### This is the Main Function and contains a simulation case

### Also CHECK THE TIME REQUIRED FOR THE MODEL
rm(list = ls())
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
library(lattice)
library(MixSim)
library(mclust)
#################################### SIMULATED DATA PROPERTIES ####################################################
options(error = recover)

source("sourceAll.R")
path = getwd()
sourceDir(path, except=c("SimData.R", "Compare.R"))
source('simulate.R')

Sim_GHDP <- function(i, tau, k){
## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0

###### Get the Data #####################################
## Initialize the Training Data
s <- 1
t <- 2
w <- 50
m <- k
p <- array(c(250),c(s,t))
n <- 70*s
r <- rep(n/s,s)

set.seed(i)
true_parm <- gen.clust(s, t, w, m, p, n, r)
X.sim <- simulate(s, t, w, m, p, n, r, tau, prob.overlap, true_parm)

### RUN
n.burn = 10000
n.reps = 50000
max.row.nbhd.size = round(.1*25*250^.5) # should be small compared to n2*p^d (~ n2*G if d=.5)
max.col.nbhd.size = round(.05*250) # should be small compared to p
row.frac.probes = 0.05
col.frac.probes = .1
prob.compute.col.nbhd=.2
prop.X.miss <- 0
tBB_flag=FALSE

data1 <- NULL
data1$X[[1]] <- X.sim[[1]][[1]]
data1$X[[2]] <- X.sim[[1]][[2]]
data <- data1

tmp <- list()

tmp[[1]] <- 1:cumsum(r)[1]
if(s >1){
  for(rr in 1:(m-1)){
    tmp[[rr+1]] <- (cumsum(r)[rr]+1):cumsum(r)[rr+1]
  }
}

All.Stuff <- fn.mcmc(s, t, m, k, p, n, r, true, data, data1, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, 
                     row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, tBB_flag, prop.X.miss)

All.Stuff

}

tau.l <- c(1, 2, 4, 5)
k.l <- c(3, 4, 5)

for(tau in tau.l){
  for (k in k.l) {
    for(i in 1:100){
      GHDP <- Sim_GHDP(i, tau, k)
      save.image(file = paste("GHDP_0_k", k, "tau", tau, "i", i, ".RData", sep = "_"))
      print(paste("GHDP_0_k", k, "tau", tau, "i", i, "DONE!", sep = "_"))
    }
  }
}


