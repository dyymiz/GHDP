
fn1.update.element.objects <- function(parm)
	{
  
	parm$clust$s.mt <- array(parm$clust$s.v, c(parm$n2, parm$clust$G))

	for (g in 1:parm$clust$G)
		{s.g.v <- parm$clust$s.mt[,g]
		s.pos.indx <- s.g.v > 0
		#
		if (sum(s.pos.indx) > 0)
			{parm$clust$A.mt[s.pos.indx,g] <- parm$clust$phi.v[s.g.v[s.pos.indx]]
			}
		if ((parm$n2-sum(s.pos.indx)) > 0)
			{parm$clust$A.mt[!s.pos.indx,g] <- 0
			}
		# parm$clust$B.mt[,(g+1)] <- parm$clust$A.mt[,g]
		}

	parm$clust$theta.v <- as.vector(parm$clust$A.mt)


	  parm$clust$n.vec <- array(,parm$clust$K)
	  for (s in 1:parm$clust$K)
		  {parm$clust$n.vec[s] <- sum(parm$clust$s.v == s)
	  }
	  parm$clust$n0 <- sum(parm$clust$s.v == 0)



	parm
	}


fn2.update.element.objects <- function(parm)
	{

  # parm$Y/S.xd is R by G, parm$n2 = parm$clust$R
	parm$Y <- parm$X.sd <- array(,c(parm$n2, parm$clust$G))

	# group covariate tells which parm$clust$rho.g
	# to use for likelihood calculation
	parm$g <- rep(1:parm$clust$G,each = parm$n2)

	for (g in 1:parm$clust$G)
		{I.g <- (parm$clust$c.v==g)
		for (r in 1:parm$clust$R) {
		  J.g <- (parm$clust$r.v==r)
		  m.g <- parm$clust$C.m.vec[g]
		  n.g <- parm$clust$R.m.vec[r]

		x.g.v <- x.tmp <- parm$X[J.g,I.g]
		x2.g.v <- x.g.v^2
		if (m.g*n.g > 1){
			x.g.v <- mean(x.tmp)
			x2.g.v <- mean(x.tmp^2)
		}
		 parm$Y[r,g] <- x.g.v

		sd.g.v <- 0
		 if ((m.g*n.g) > 1){
			sd.g.v <- sqrt((x2.g.v - x.g.v^2)*(m.g*n.g)/((m.g*n.g)-1))
		 }
		parm$X.sd[r,g] <- sd.g.v
	}
		}

	parm$clust$R_C <- c(outer(parm$clust$R.m.vec , parm$clust$C.m.vec))
	
	parm$N <- parm$n2 * parm$clust$G

	parm$Y <- as.vector(parm$Y)

	parm$X.sd <- as.vector(parm$X.sd)

	#####################

	parm <- fn1.update.element.objects(parm)

	parm
	}



fn.element.DP <- function(data, parm, max.row.nbhd.size, row.frac.probes)
{

	# essentially, a Bush-Mac move: given groups, the parm$N=n2XG number of
	# invidividual elements (summaries of microarray elements) belonging to group g>0
	# are updated for s (phi) and z

	parm <- fn2.update.element.objects(parm)

	parm <- element_fn.fast.DP(parm, max.row.nbhd.size, row.frac.probes)

	#############################
	## Important: do not remove call to fn1.update.element.objects
	## updates A.mt, theta.v, B.mt, tBB.mt, s.v, s.mt, n.vec, n0
	#############################

	parm <- fn1.update.element.objects(parm)

  	parm
}



