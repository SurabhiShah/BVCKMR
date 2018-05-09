########################
# Data simulation code #
########################

library(MASS); library(stats)

#n 				= 100 # Number of subjects
q 				= 2 # Random slope and random intercept
p 				= 2 # Number of confounders
M				  = 5 # Number of metals
age 			= c(-2:2)
T				  = length(age) # Number of time points
N				  = n*T

#############################

# Simulating Z as multivariate normal matrix

Z 				= matrix(mvrnorm(n = n, mu=c(rep(0,M)), Sigma = matrix(data = c(1,0.2,0.3,0.45,0.15, 0.2,1,0.4,0.45,0.3, 0.3,0.4,1,0.3,0.15, 0.45,0.45,0.3,1,0.2, 0.15,0.3,0.15,0.2,1), nrow=M)), byrow=FALSE, nrow=n)

Z 				= apply(Z, 2, scale)

X				= cbind( matrix( rep(rnorm(n, mean=0, sd=1), each=T), byrow=FALSE ), matrix( rep(scale(matrix(sample(1:2, n, replace = TRUE))), each=T) ))
beta			= rep(1,p)

b				= c( t(mvrnorm(n = n, mu=rep(0,q), Sigma = matrix(data = c(0.25,0.005,0.005,0.05), nrow=q))))

res.sd = 1

# Make U matrix

U				= matrix(data=rep(0,T*2*n), nrow=T)
U[,1]			= rep(1, T)
U[,2]			= age

counter 		= 1
while(counter < n) {
    U_i					= matrix(data=rep(0,T*2*n), nrow=T)
    U_i[,counter*2+1]	= rep(1, T)
    U_i[,counter*2+2]	= age
    U 					= rbind(U, U_i)
    counter 				= counter+1
}

Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep(0.5*(Z[,1]^2 + Z[,2]^2 - Z[,4]^2 - Z[,5]^2 + 0.5*Z[,4]*Z[,5] + Z[,4] + Z[,5]), each=T) + matrix(rep(0.5*(Z[,1]^2 - Z[,2]^2 - Z[,1]*Z[,2] + Z[,3]^2 + Z[,4] - Z[,5]), each=T)) * rep(age, n) + U %*% b

# Make W matrix

W				= matrix(data=rep(0,T*2*n), nrow=T)
W[,1]			= rep(1, T)
W[,n+1]			= age

counter 		= 2
while(counter <= n) {
	W_i				= matrix(data=rep(0,T*2*n), nrow=T)
	W_i[,counter]	= rep(1, T)
	W_i[,n+counter]	= age
	W 				= rbind(W, W_i)
	counter 		= counter+1
}

WTW = t(W) %*% W
XTX	= t(X) %*% X
UTU = t(U) %*% U
