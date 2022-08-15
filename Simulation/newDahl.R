
fn.Dahl <- function(parm, All.Stuff, dahl.flag, mm)
{### Dahl calculations
  n2 <- dim(parm$X)[2]
  tmp.mat <- array(0,c(n2,n2))
  
  for (jj in 1:parm$clust$G)
  {indx.jj <- which(parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  parm$dahlDist.mt = tmp.mat
  
  if (dahl.flag)
  {# interpretation: square of this is double the total number of mismatching pairs of allocations
    parm$Dahl_dist =  sqrt(sum((parm$dahlDist.mt - All.Stuff$meanDahlDist.mt[[mm]])^2))
    
  if (parm$Dahl_dist < parm$min_Dahl_dist)
  {parm$min_Dahl_dist = parm$Dahl_dist
  parm$est_c.v  <- parm$clust$c.v
  }
  }
  
  parm
}



fn.Dahl.c <- function(parm, All.Stuff, dahl.flag, mm, tt)
{### Dahl calculations

  tmp.mat <- array(0,c(parm$p,parm$p))
  
  for (jj in 1:parm$clust$G)
  {indx.jj <- which(parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  parm$dahlDist.mt = tmp.mat
  
  if (dahl.flag)
  {# interpretation: square of this is double the total number of mismatching pairs of allocations
    parm$Dahl_dist =  sqrt(sum((parm$dahlDist.mt - All.Stuff$meanDahlDist.mt.c[[mm]][[tt]])^2))
    
    if (parm$Dahl_dist < parm$min_Dahl_dist)
    {parm$min_Dahl_dist = parm$Dahl_dist
    parm$est_c.v  <- parm$clust$c.v
    }
  }
  
  parm
}
