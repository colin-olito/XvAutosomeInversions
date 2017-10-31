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
source('R/functions-RecSim-Autosomal.R')


par.list  <-  list(
					gen  =  25000,
					m    =  0.01,
					s    =  0.1,
					h    =  0,
					r    =  0.05
				   )

xi.init  <-  c(0,(par.list$m/par.list$s),(par.list$m/par.list$s),(1 - 2*(par.list$m/par.list$s)-0.01),0.01)
res  <-  recursionFwdSim(par.list, xi.init = xi.init, threshold = 1e-6)

str(res)

head(res$xi.gen)

plot(NA, ylim=c(0,1), xlim=c(0,nrow(res$xi.gen)), ylab="Frequency", xlab="Generations")
	lines(res$xi.gen[,1], col=1, lwd=3)
	lines(res$xi.gen[,2], col=8, lwd=3, lty=1)
	lines(res$xi.gen[,3], col=4, lwd=3, lty=3)
	lines(res$xi.gen[,4], col=3, lwd=3)
	lines(res$xi.gen[,5], col=2, lwd=3)
res$EQ.freq