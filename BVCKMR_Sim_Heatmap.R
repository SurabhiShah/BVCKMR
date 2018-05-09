###########################
# Heatmap simulation plot #
###########################

source("BVCKMR_Sim_Pred.R")

ngrid <- 30
qs <- c(0.25, 0.5, 0.75, 0.05, 0.95)
j = 2 # median levels 
min.plot.dist <- 0.5 ## only color level plot for points within this distance from an observed data point

ylim = c(-5, 5)

######################################################
# Estimated and true exposure-response relationships #
######################################################

z1.1 <- seq(quantile(Z[,1], qs[4]), quantile(Z[,1], qs[5]), length=ngrid)
z1.2 <- seq(quantile(Z[,2], qs[4]), quantile(Z[,2], qs[5]), length=ngrid)

cross.sec.1 = cbind(expand.grid(z1=z1.1, z2 = z1.2), z3=median(Z[,3]), z4=median(Z[,4]), z5=median(Z[,5]))

grid = expand.grid(z1=z1.1, z2 = z1.2)
mindists = rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
		pt <- as.numeric(grid[i,])
		dists <- as.matrix(dist(rbind(pt, Z[,1:2])))["pt",-1]
		mindists[i] <- min(dists)
}
rm(grid, pt, dists)

# Specific cross-sectional graphs: 

# Baseline, h1
hgrid.T.1 <- newh.postmean(Znew=cross.sec.1, sel=sel)$postmean[(ngrid*ngrid+1):(2*ngrid*ngrid)]
hgrid.T.1[mindists > min.plot.dist] <- NA
hgrid.T.1 = matrix(hgrid.T.1, nrow=ngrid)
# Trajectory, h2
hgrid.T.2 = newh.postmean(Znew=cross.sec.1, sel=sel)$postmean[1:(ngrid*ngrid)]
hgrid.T.2[mindists > min.plot.dist] <- NA
hgrid.T.2 = matrix(hgrid.T.2, nrow=ngrid)

# Truth:
  true.1 = matrix(0.5*(cross.sec.1[,1]^2 + cross.sec.1[,2]^2 - cross.sec.1[,4]^2 - cross.sec.1[,5]^2 + 0.5*cross.sec.1[,4]*cross.sec.1[,5] + cross.sec.1[,4] + cross.sec.1[,5]), nrow=ngrid)
  true.2 = matrix(0.5*(cross.sec.1[,1]^2 - cross.sec.1[,2]^2 - cross.sec.1[,1]*cross.sec.1[,2] + cross.sec.1[,3]^2 + cross.sec.1[,4] - cross.sec.1[,5]), nrow=ngrid)

######################################################################################
# Plotting a panel of the estimated LKMR heatmap (top) and the true heatmap (bottom) #
######################################################################################

#### multiple images with a common legend
#pdf("Plot.pdf", height=6, width=6)

set.panel()
zlim = c(-2,4)

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) # 2X2 matrix of plots

# now draw all your plots using usual image command

image(z1.1, z1.2, hgrid.T.1, xlab="z1", ylab="z2", col=tim.colors(), main="BVCKMR Estimated h1(z1, z2)", zlim=zlim)
#points(Z)
image(z1.1, z1.2, hgrid.T.2, xlab="z1", ylab="z2", col=tim.colors(), main="BVCKMR Estimated h2(z1, z2)", zlim=zlim)
#points(Z)

image(z1.1, z1.2, true.1, xlab="z1", ylab="z2", col=tim.colors(), main="True h1(z1, z2)", zlim=zlim)
#points(Z)
image(z1.1, z1.2, true.2, xlab="z1", ylab="z2", col=tim.colors(), main="True h2(z1, z2)", zlim=zlim)
#points(Z)

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=zlim) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device
#dev.off()

