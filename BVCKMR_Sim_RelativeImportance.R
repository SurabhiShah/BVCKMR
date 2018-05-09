############################################################
# Relative importance at baseline (h1) and trajectory (h2) #
############################################################

source("BVCKMR_Sim_Pred.R")

par(mfrow=c(1,2))

mat 		= matrix(NA, nrow=q, ncol=M)
mat.sd 		= matrix(NA, nrow=q, ncol=M)
true.mat	= matrix(rep(0,q*M), nrow=q, ncol=M)

qs = c(0.25, 0.75, 0.50)

for (j in 1:M) {
		cross.sec 		= rbind(apply(Z, 2, median), apply(Z, 2, median))
		cross.sec[,j] 	= c(quantile(Z[,j], qs[2]), quantile(Z[,j], qs[1]))
		hgrid.cross.sec = newh.postmean(Znew = cross.sec, sel = sel) #returns 75% h2, 25% h2, 75% h1, 25% h1
		mat[1,j] 	= hgrid.cross.sec$postmean[3] - hgrid.cross.sec$postmean[4]
		mat[2,j] 	= hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
		mat.sd[1,j] 	= sqrt(hgrid.cross.sec$postvar[3,3] + hgrid.cross.sec$postvar[4,4] - 2*hgrid.cross.sec$postvar[3,4])
		mat.sd[2,j] 	= sqrt(hgrid.cross.sec$postvar[1,1] + hgrid.cross.sec$postvar[2,2] - 2*hgrid.cross.sec$postvar[1,2])
	}
	
#Plot
ylim 		= c(-1, 1)
plot(1:M, mat[1,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Baseline (h1)"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(1:M, mat[1,] - 1.96*mat.sd[1,], 1:M, mat[1,] + 1.96*mat.sd[1,], length=0.05, angle=90, code=3)
axis(1, at=1:M, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)

#Plot
ylim 		= c(-1, 1)
plot(1:M, mat[2,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Trajectory (h2)"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(1:M, mat[2,] - 1.96*mat.sd[2,], 1:M, mat[2,] + 1.96*mat.sd[2,], length=0.05, angle=90, code=3)
axis(1, at=1:M, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)
