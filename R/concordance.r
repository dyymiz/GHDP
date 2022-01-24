

num.indx <- length(indx.set)
num.indx2 <- length(indx2.set)

count <- num.permissible <- num.concordant <- 0

for (i in 1:num.indx)
for (j in 1:num.indx2)
	{
	i2 <- sort(indx.set)[i]
	j2 <- sort(indx2.set)[j]
	k <- ifelse((data$true$Y[i2] < data$true$Y[j2]), i2, j2)
	K <- i2 + j2 - k	
	pred.k <- mean.v[k] 
	pred.K <- mean.v[K] 

	if ((data$true$delta[k] == 1)&(k != K))
		{count <- count + 1

		num.permissible <- num.permissible + 1
		if (pred.k < pred.K)
			{num.concordant <- num.concordant + 1
			}
		if (pred.k == pred.K)
			{num.concordant <- num.concordant + .5
			}
		}
	}

concord.indx <- num.concordant/num.permissible
error.rate <- 1 - concord.indx

############################


