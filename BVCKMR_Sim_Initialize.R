##################################
# Simulation initialization code #
##################################

library(lars); library(lasso2); library(mvtnorm); library(SuppDists); library(MCMCpack); 
library(grplasso); library(magic); library(kernlab); library(MASS); library(fields); library(stats)

# Set notation. 
# T = number of follow-up time points
# n = Number of subjects
# N = T*n, assuming that every subject is followed up the same number of instances
# q = 2 for random slope and random intercept
# p = number of confounders
# Z = metal matrix
# U = random effects matrix

sig.shape 	= 85 
sig.scale 	= 10
sig.sq		  = rinvgamma(1, shape = sig.shape, scale = sig.scale)

v0			    = q
C0			    = diag(q)
D.inv		    = matrix(rWishart(1, df = v0, Sigma=C0), nrow=q)

r1			    = 1
delta1		  = 1
lambda1.sq 	= rgamma(1,shape = r1, rate = delta1) 
tau.sq  	  = rgamma(q, shape = (n + 1)/2, rate = lambda1.sq/2)

beta		    = rep(1, p)
beta.num	  = length(beta)

poly        = polydot(degree=2, offset=1)

# Sigma_h matrix
list.K				  = list()
poly 				    = polydot(degree=2, offset=1)
for (k in 1:q) {
	list.K[[k]] 	= kernelMatrix(poly, Z)*tau.sq[k]
	}  
Sigma.h 			    = do.call(adiag, list.K)
if (doginv) {
	Sigma.h.inv		  = ginv(Sigma.h)
	} else if (docofactor) {
		Sigma.h.inv		= solve(Sigma.h + diag(q*n) * cofactor)
		} else Sigma.h.inv 	= solve(Sigma.h) 
mean.h 				    = rep(0, q*n)    			

h  					      = mvrnorm(1, mu=mean.h, Sigma=Sigma.h)  
    
# FOR POSTERIOR
sigsq.post 	= lambda1.post = NULL 
D.inv.post 	= list()
b.post		  = rbind( b, matrix(rep(NA, num.reps*2*n), ncol=2*n) )
h.post      = rbind( h, matrix(rep(NA,num.reps*2*n), ncol=2*n) )
tausq.post 	= rbind( tau.sq, matrix(rep(NA,num.reps*q), ncol=q) )
beta.post 	= rbind( beta, matrix(rep(NA,num.reps*beta.num), nrow=num.reps) ) 
h			      = matrix(h, nrow=2*n)



