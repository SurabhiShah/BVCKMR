##############################
# BVCKMR Simulation Run Code #
##############################

# Set tasks
doMCMC 			= TRUE
doginv			= FALSE
docofactor	= TRUE
cofactor		= 1e-7	# randomly generated kernel matrix is not always invertible.

if (doMCMC)
{
	n			    = 100
	source("BVCKMR_Sim_Data.R")
  num.reps 	= 20000
	sel 		  = 10000:num.reps
	
	#cat(c("Job started at:",date()),fill=TRUE)
	source("BVCKMR_Sim_Initialize.R")
	source("BVCKMR_MCMC.R")
	#cat(c("Job finished at:",date()),fill=TRUE)
	MCMC = list(Sigsq = sigsq.post, Lam1 = lambda1.post, H = h.post, Tau = tausq.post, Beta = beta.post, DINV = D.inv.post, b = b.post, Z = Z, Y = Y, X = X, U = U, W = W, sel = sel, cofactor = cofactor, q = q, M = M, n = n)  
	
	#save(MCMC, file=paste0("MCMC", jobid, ".RData"))
}

# Analysis

Summary.h 	= matrix(NA, nrow=q, ncol=4)

true.h1 = 0.5*(Z[,1]^2 + Z[,2]^2 - Z[,4]^2 - Z[,5]^2 + 0.5*Z[,4]*Z[,5] + Z[,4] + Z[,5])
true.h2 = 0.5*(Z[,1]^2 - Z[,2]^2 - Z[,1]*Z[,2] + Z[,3]^2 + Z[,4] - Z[,5])

hhat.1	= apply(MCMC$H[sel, 1:n], 2, mean)
hhat.2	= apply(MCMC$H[sel, (n+1):(2*n)], 2, mean)

model.1 = lm(hhat.1~true.h1)
model.2 = lm(hhat.2~true.h2)

Summary.h[1,] = c(summary(model.1)$coef[1,1], summary(model.1)$coef[2,1], summary(model.1)$r.sq, sqrt( sum( (hhat.1 - true.h1)^2) / n))
Summary.h[2,] = c(summary(model.2)$coef[1,1], summary(model.2)$coef[2,1], summary(model.2)$r.sq, sqrt( sum( (hhat.2 - true.h2)^2) / n))

