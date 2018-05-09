for (iter in 1:num.reps)  {
	# Sigma_h matrix
	list.K			= list()
	poly 			= polydot(degree=2, offset=1)
	for (k in 1:q) {
		list.K[[k]] = kernelMatrix(poly, Z)*tau.sq[k]
	}  
	Sigma.h 		= do.call(adiag, list.K)
	if (doginv) {
		Sigma.h.inv		= ginv(Sigma.h)
		} else if (docofactor) {
			Sigma.h.inv		= solve(Sigma.h + diag(q*n) * cofactor)
			} else Sigma.h.inv 	= solve(Sigma.h) 
	cov.h 			= solve( WTW / sig.sq + Sigma.h.inv)
    mean.h    		= solve( WTW / sig.sq + Sigma.h.inv) %*% t(W) %*% (Y - X %*% beta - U %*% b) * (1/sig.sq)
    h				= mvrnorm(1,mu=mean.h,Sigma=cov.h)
    h.post[iter+1,] = h
    h 				= matrix(h, nrow=2*n)

    # sig.sq
    sh.sig      	= N/2 + n + sig.shape 
    kron.D.inv		= diag(n) %x% D.inv
    sc.sig      	= 1/2*t(Y - X %*% beta - W %*% h - U %*% b) %*% (Y - X %*% beta - W %*% h - U %*% b) + 1/2*t(b) %*% kron.D.inv %*% b + sig.scale
    sig.sq     		= rinvgamma(1, shape=sh.sig, scale=sc.sig)
    sigsq.post 		= c(sigsq.post, sig.sq)
     
    gam 			= c()
    for (j in 1:q){
    	term.h 		= h[((j-1)*n+1):(j*n)]
    	if (doginv) {
		term.gam		= t(term.h) %*% ginv(kernelMatrix(poly, Z)) %*% term.h
			} else if (docofactor) {
			term.gam= t(term.h) %*% solve(kernelMatrix(poly, Z) + 							diag(n)*cofactor) %*% term.h
			} else term.gam = t(term.h) %*% solve(kernelMatrix(poly, Z)) %*% 					term.h
    	gam[j]  	= rinvGauss(n=1, nu=sqrt(lambda1.sq*sig.sq/term.gam), lambda=lambda1.sq)            		
    	tau.sq[j] 	= 1/gam[j] 
    }  	
    tausq.post[iter+1,] = tau.sq 
	    
    sh.lam1       = n + 1 + r1 
    sc.lam1       = 1/2*sum(tau.sq) + delta1
    lambda1.sq    = rgamma(1, shape=sh.lam1, rate=sc.lam1)
    lambda1.post  = c(lambda1.post, lambda1.sq)
    
    mean.beta 			= solve(XTX) %*% t(X) %*% (Y - W %*% h - U %*% b)
    sig.beta 			= sig.sq * solve(XTX) 
    beta 				= mvrnorm(1, mu=mean.beta, Sigma=sig.beta)
    beta.post[iter+1,] 	= beta 
    beta 				= matrix(beta, nrow=p)
    
    b.hat				= solve(UTU) %*% t(U) %*% (Y - X %*% beta - W %*% h)
    b.star				= solve( (diag(n) %x% D.inv) + UTU ) %*% (UTU %*% b.hat)
    mean.b				= b.star
    sig.b				= sig.sq * solve( (diag(n) %x% D.inv) + UTU )
    b					= mvrnorm(1, mu=mean.b, Sigma=sig.b)
    b.post[iter+1,]		= b
    b 					= matrix(b, nrow=2*n)
    
    df.D.inv			= n + v0
    B.D					= cbind(b[c(TRUE, FALSE)], b[c(FALSE, TRUE)])
    Sigma.D.inv			= solve( solve(C0) + t(B.D) %*% B.D / sig.sq)
    D.inv				= matrix(rWishart(1, df=df.D.inv, Sigma = Sigma.D.inv), nrow=q)
    D.inv.post[[iter]]	= D.inv
}
