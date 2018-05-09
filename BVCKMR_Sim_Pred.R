#####################################
# Prediction code for new exposures #
# Postmean estimation functions     #
#####################################

newh.postmean <- function(Znew, sel) {
	
	Z = MCMC$Z
	
	if(is.null(dim(Znew))) Znew = matrix(Znew, nrow=1)
	if(class(Znew)  == "data.frame") Znew <- data.matrix(Znew)
	Znew            = as.matrix(Znew)
	pnew            = nrow(Znew)
	
	Sigsq.mean	    = mean(MCMC$Sigsq[sel])
	Tau.mean	      = apply(MCMC$Tau[sel,], 2, mean)
	b.mean		      = apply(MCMC$b[sel,], 2, mean)
	Beta.mean	      = apply(MCMC$Beta[sel,], 2, mean)
	h.mean		      = apply(MCMC$H[sel,], 2, mean)
	Y			          = MCMC$Y
	X 			        = MCMC$X
	U			          = MCMC$U
	W               = MCMC$W
		
	list.K				  = list()

		poly 				    = polydot(degree=2, offset=1)
	for (k in 1:q) {
		list.K[[k]] 	= kernelMatrix(poly, Z)*Tau.mean[k]
	}  
	Sigma.h 				= do.call(adiag, list.K)
	
	Sigma.h.new = "[<-"(matrix(0, (q*n+q*pnew), (q*n+q*pnew)), 1:nrow(Sigma.h), 1:ncol(Sigma.h), value = Sigma.h) # adding rows/cols of zeroes to square matrix

  list.G        = list()
  for (ind in 1:q) {
		list.G[[ind]] = kernelMatrix(poly, rbind(Z, Znew))*Tau.mean[ind]
	}   

	# h vector is ordered c(h1, h2, h2new, h1new)	
	G.1.11 = list.G[[1]][ 1:n, 1:n ]
	G.1.22 = list.G[[1]][ (n+1):(n+pnew), (n+1):(n+pnew) ]
	G.1.21 = list.G[[1]][ (n+1):(n+pnew), 1:n ]
	G.1.12 = list.G[[1]][ 1:n, (n+1):(n+pnew) ]
	
	list.G.final = list()
	list.G.final[[1]] = G.1.11
	list.G.final[[2]] = list.G[[2]]
	list.G.final[[3]] = G.1.22	
	cov.bf = do.call(adiag, list.G.final)
	
	cov.bf[(q*n + pnew + 1):(q*n + q*pnew), 1:n] = G.1.21
	cov.bf[1:n, (q*n + pnew + 1):(q*n + q*pnew)] = G.1.12
	
	Sigma11 = cov.bf[1:(n*q), 1:(n*q)]
	Sigma12 = cov.bf[1:(n*q), (n*q + 1):(n*q + pnew*q)] 
	Sigma21 = cov.bf[(n*q + 1):(n*q + pnew*q), 1:(n*q)]
	Sigma22 = cov.bf[(n*q + 1):(n*q + pnew*q), (n*q + 1):(n*q + pnew*q)]
	
	WTW         = t(W) %*% W
	Sigma.h.inv = solve(Sigma.h + cofactor*diag(nrow(Sigma.h)))
	m.h         = solve( WTW / Sigsq.mean + Sigma.h.inv) %*% t(W) %*% (Y - X %*% Beta.mean - U %*% b.mean) * (1/Sigsq.mean)
	s.h         = solve( WTW / Sigsq.mean + Sigma.h.inv)
	
	Sigma11.inv = Sigma.h.inv 
	
	Mu.hnew.updated = Sigma21 %*% Sigma11.inv %*% m.h
	Sigma.hnew.updated = Sigma22 - Sigma21 %*% Sigma11.inv %*% Sigma12 + Sigma21 %*% Sigma11.inv %*% s.h %*% Sigma11.inv %*% Sigma12 

	list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))
}
