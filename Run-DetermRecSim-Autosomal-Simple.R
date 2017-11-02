##############################################################
#  Local adaptation and the evolution of autosomal inversions
#
#  R code for simple deterministic simulations of the
#  haplotype frequency recursions for the model of 
#  autosomal inversions. Simplest case to find equilibrium
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
source('R/functions-DetermRecSim-Autosomal.R')


######################
##  Run Simulations

# Locally adaptive allele completely recessive (h = 0)
test  <-  recursionFwdSimLoop(gen = 25000, h = 0, s.vals = c(0.05, 0.1, 0.2),
                    m.vals = c(0.01, 0.05), r.vals = c(0.0, 0.01), 
                    threshold = 1e-7)

# Additive fitness effects (h = 1/2)
test  <-  recursionFwdSimLoop(gen = 25000, h = 0.5, s.vals = c(0.05, 0.1, 0.2),
                    m.vals = c(0.01, 0.05), r.vals = c(0.0, 0.01), 
                    threshold = 1e-7)

# Locally adaptive allele completely dominant (h = 1)
test  <-  recursionFwdSimLoop(gen = 25000, h = 1, s.vals = c(0.05, 0.1, 0.2),
                    m.vals = c(0.01, 0.05), r.vals = c(0.0, 0.01), 
                    threshold = 1e-7)


########################################################
##  Some exploratory code to play with results/plotting

head(test)
unique(test$s)
rs  <-  unique(test$r)
ms  <-  unique(test$m)

 plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(test$s)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
       	lines(test$x5[test$r==rs[1] & test$m==ms[1]] ~ test$s[test$r==rs[1] & test$m==ms[1]], col=1, lwd=3, lty=1)
       	lines(test$x5[test$r==rs[1] & test$m==ms[2]] ~ test$s[test$r==rs[1] & test$m==ms[2]], col=1, lwd=3, lty=1)
       	lines(test$x5[test$r==rs[2] & test$m==ms[1]] ~ test$s[test$r==rs[2] & test$m==ms[1]], col=2, lwd=3, lty=2)
       	lines(test$x5[test$r==rs[2] & test$m==ms[2]] ~ test$s[test$r==rs[2] & test$m==ms[2]], col=2, lwd=3, lty=2)
 # axes
        axis(1, las=1)
        axis(2, las=1)


 plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(test$s)), ylim = c(0,max(test$LD)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(test$LD[test$r==rs[1] & test$m==ms[1]] ~ test$s[test$r==rs[1] & test$m==ms[1]], col=1, lwd=3, lty=1)
        lines(test$LD[test$r==rs[1] & test$m==ms[2]] ~ test$s[test$r==rs[1] & test$m==ms[2]], col=2, lwd=3, lty=1)
        lines(test$LD[test$r==rs[2] & test$m==ms[1]] ~ test$s[test$r==rs[2] & test$m==ms[1]], col=1, lwd=3, lty=2)
        lines(test$LD[test$r==rs[2] & test$m==ms[2]] ~ test$s[test$r==rs[2] & test$m==ms[2]], col=2, lwd=3, lty=2)
 # axes
        axis(1, las=1)
        axis(2, las=1)


########################################################
##  Some exploratory code to play with results/plotting
#
# par.list  <-  list(
# 					gen  =  25000,
# 					m    =  0.01,
# 					s    =  0.1,
# 					h    =  0.5,
# 					r    =  0.05
# 				   )
#
# xi.init  <-  c(0,(par.list$m/par.list$s),(par.list$m/par.list$s),(1 - 2*(par.list$m/par.list$s)-0.01),0.01)
## xi.init  <-  c(0.25,0.25,0.25,(0.25-0.01),0.01)
#
# res  <-  recursionFwdSim(par.list, xi.init = xi.init, threshold = 1e-6, verbose=FALSE)
# str(res)
# head(res$xi.gen)
#
# plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$xi.gen)), ylab="Frequency", xlab="Generations")
# 	lines(res$xi.gen[,1], col=1, lwd=3)
# 	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
# 	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
# 	lines(res$xi.gen[,4], col=3, lwd=3)
# 	lines(res$xi.gen[,5], col=2, lwd=3)
# res$EQ.freq
#
#
## Finding equilibrium haplotype frequencies in the absence of the inversion,
## then using these eq. frequencies as initial conditions when inversion invades
#
# inits  <-  recursionFwdSim(par.list, xi.init = c(0.25,0.25,0.25,0.25,0), threshold = 1e-6)
# inits  <-  inits[[3]]
# inits[inits == max(inits)]  <-  inits[inits == max(inits)] - 0.01
# inits[5]  <-  0.01
#
# res  <-  recursionFwdSim(par.list, xi.init = inits, threshold = 1e-7)
#
# plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$xi.gen)), ylab="Frequency", xlab="Generations")
# 	lines(res$xi.gen[,1], col=1, lwd=3)
# 	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
# 	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
# 	lines(res$xi.gen[,4], col=3, lwd=3)
# 	lines(res$xi.gen[,5], col=2, lwd=3)
# res$EQ.freq
