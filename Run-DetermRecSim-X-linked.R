##############################################################
#  Local adaptation and the evolution of X-linked inversions
#  with sex-specific selection and migration
#
#  R code for simple deterministic simulations of the
#  haplotype frequency recursions for the model of 
#  X-linked inversions with sex-specific selection and 
#  migration. Simplest case to find equilibrium
#  frequencies of the inversions.
#
#
#  Author: Colin Olito
#
#  NOTES:  
#          

rm(list=ls())
#####################
##  Dependencies

source('R/functions-figures.R')
source('R/functions-DetermRecSim-X-linked.R')



######################
##  Run Simulations

# Locally adaptive allele completely recessive (hf = hm = 0)
test  <-  recursionFwdSimLoop(gen = 25000, h = 0, threshold = 1e-7,
					sf.vals = c(0.05, 0.1, 0.2), sm.vals = c(0.05, 0.1, 0.2),
					mf.vals = c(0.01, 0.05), mm.vals = c(0.01, 0.05),
					r.vals = c(0.0, 0.01, 0.1))

# Additive fitness effects  hf = hm = 1/2)
test  <-  recursionFwdSimLoop(gen = 25000, h = 0.5, threshold = 1e-7,
					sf.vals = c(0.05, 0.15, 0.2), sm.vals = c(0.05, 0.1, 0.2),
					mf.vals = c(0.01, 0.05), mm.vals = c(0.01, 0.05),
					r.vals = c(0.0, 0.01, 0.1))

# Locally adaptive allele completely dominant (hf = hm = 1)
test  <-  recursionFwdSimLoop(gen = 25000, h = 1, threshold = 1e-7,
					sf.vals = c(0.01, 0.05, 0.1, 0.2), sm.vals = c(0.01, 0.05, 0.1, 0.2),
					mf.vals = c(0.01, 0.05), mm.vals = c(0.01, 0.05),
					r.vals = c(0.0, 0.01, 0.1))



########################################################
##  Some exploratory code to play with results/plotting
head(test)

unique(test$sf)
rs       <-  unique(test$r)
mfs      <-  unique(test$mf)
EQ.freq  <-  (2*test[,8:12] + test[,13:17])/2

 plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(test$sf)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
       	lines(test$x5[test$r==rs[1] & test$mf==mfs[1]] ~ test$sf[test$r==rs[1] & test$mf==mfs[1]], col=1, lwd=3, lty=1)
       	lines(test$x5[test$r==rs[1] & test$mf==mfs[2]] ~ test$sf[test$r==rs[1] & test$mf==mfs[2]], col=1, lwd=3, lty=1)
       	lines(test$x5[test$r==rs[2] & test$mf==mfs[1]] ~ test$sf[test$r==rs[2] & test$mf==mfs[1]], col=1, lwd=3, lty=1)
       	lines(test$x5[test$r==rs[2] & test$mf==mfs[2]] ~ test$sf[test$r==rs[2] & test$mf==mfs[2]], col=1, lwd=3, lty=1)
       	lines(test$y5[test$r==rs[1] & test$mf==mfs[1]] ~ test$sf[test$r==rs[1] & test$mf==mfs[1]], col=2, lwd=3, lty=2)
       	lines(test$y5[test$r==rs[1] & test$mf==mfs[2]] ~ test$sf[test$r==rs[1] & test$mf==mfs[2]], col=2, lwd=3, lty=2)
       	lines(test$y5[test$r==rs[2] & test$mf==mfs[1]] ~ test$sf[test$r==rs[2] & test$mf==mfs[1]], col=2, lwd=3, lty=2)
       	lines(test$y5[test$r==rs[2] & test$mf==mfs[2]] ~ test$sf[test$r==rs[2] & test$mf==mfs[2]], col=2, lwd=3, lty=2)
 # axes
        axis(1, las=1)
        axis(2, las=1)

 plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(test$sf)), ylim = c(0,max(test$LDx)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
       	lines(test$LDx[test$r==rs[1] & test$mf==mfs[1]] ~ test$sf[test$r==rs[1] & test$mf==mfs[1]], col=1, lwd=3, lty=1)
       	lines(test$LDx[test$r==rs[1] & test$mf==mfs[2]] ~ test$sf[test$r==rs[1] & test$mf==mfs[2]], col=1, lwd=3, lty=1)
       	lines(test$LDx[test$r==rs[2] & test$mf==mfs[1]] ~ test$sf[test$r==rs[2] & test$mf==mfs[1]], col=1, lwd=3, lty=1)
       	lines(test$LDx[test$r==rs[2] & test$mf==mfs[2]] ~ test$sf[test$r==rs[2] & test$mf==mfs[2]], col=1, lwd=3, lty=1)
       	lines(test$LDy[test$r==rs[1] & test$mf==mfs[1]] ~ test$sf[test$r==rs[1] & test$mf==mfs[1]], col=2, lwd=3, lty=2)
       	lines(test$LDy[test$r==rs[1] & test$mf==mfs[2]] ~ test$sf[test$r==rs[1] & test$mf==mfs[2]], col=2, lwd=3, lty=2)
       	lines(test$LDy[test$r==rs[2] & test$mf==mfs[1]] ~ test$sf[test$r==rs[2] & test$mf==mfs[1]], col=2, lwd=3, lty=2)
       	lines(test$LDy[test$r==rs[2] & test$mf==mfs[2]] ~ test$sf[test$r==rs[2] & test$mf==mfs[2]], col=2, lwd=3, lty=2)
 # axes
        axis(1, las=1)
        axis(2, las=1)



#Some exploratory code to play with the recursionFwdSim function
par.list  <-  list(
			   gen  =  25000,
			   mf    =  0.05,
			   mm    =  0.01,
			   sf    =  0.2,
			   sm    =  0.1,
			   h     =  0.5,
			   r     =  0.5
			   )

xi.init  <-  c(0.25,0.25,0.25,(0.25-0.01),0.01)
yi.init  <-  c(0.25,0.25,0.25,(0.25-0.01),0.01)

res  <-  recursionFwdSim(par.list=par.list, xi.init = xi.init, yi.init = yi.init, threshold = 1e-6, silent=FALSE)

str(res)
head(res$Fi.gen)
plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$Fi.gen)), ylab="Frequency", xlab="Generations")
	lines(res$xi.gen[,1], col=1, lwd=3)
	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$xi.gen[,4], col=3, lwd=3)
	lines(res$xi.gen[,5], col=2, lwd=3)


plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$xi.gen)), ylab="Frequency", xlab="Generations")
	lines(res$xi.gen[,1], col=1, lwd=3)
	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$xi.gen[,4], col=3, lwd=3)
	lines(res$xi.gen[,5], col=2, lwd=3)
	lines(res$yi.gen[,1], col=1, lwd=3)
	lines(res$yi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$yi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$yi.gen[,4], col=3, lwd=3)
	lines(res$yi.gen[,5], col=2, lwd=3)
sum(res$EQ.freq)

res$EQ.freq
# Finding equilibrium haplotype frequencies in the absence of the inversion,
# then using these eq. frequencies as initial conditions when inversion invades

inits       <-  recursionFwdSim(par.list, xi.init = c(0.25,0.25,0.25,0.25,0), yi.init = c(0.25,0.25,0.25,0.25,0), threshold = 1e-6)
x.inits     <-  inits[[5]][1:5]
x.inits[x.inits == max(x.inits)]  <-  x.inits[x.inits == max(x.inits)] - 0.01
x.inits[5]  <-  0.01
x.inits
y.inits     <-  inits[[5]][6:10]
y.inits[y.inits == max(y.inits)]  <-  y.inits[y.inits == max(y.inits)] - 0.01
y.inits[5]  <-  0.01
y.inits
sum(x.inits)
sum(y.inits)
round(sum(y.inits), digits=3)
y.inits     <-  round(y.inits,digits=3)


res  <-  recursionFwdSim(par.list, xi.init = x.inits, yi.init = y.inits, threshold = 1e-7)
plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$xi.gen)), ylab="Frequency", xlab="Generations")
	lines(res$xi.gen[,1], col=1, lwd=3)
	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$xi.gen[,4], col=3, lwd=3)
	lines(res$xi.gen[,5], col=2, lwd=3)
	lines(res$yi.gen[,1], col=1, lwd=3)
	lines(res$yi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$yi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$yi.gen[,4], col=3, lwd=3)
	lines(res$yi.gen[,5], col=2, lwd=3)
sum(res$EQ.freq)









