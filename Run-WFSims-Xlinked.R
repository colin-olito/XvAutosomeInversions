##############################################################
#  Wright-Fisher forward simulations of invasion of an X-LINKED
#  inversion capturing 2 locally adaptive alleles, as well as 
#  possible linked deleterious mutations
#
#  R code for W-F forward simulations. Generates output data
#  as .csv files saved to ./output/data/simResults.
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		
#		About chosen parameter values:
#			Two major constraints on parameter values. 
#			1) initial freq. of maladaptive alleles should be low ((2*m)/s < 0.1) &
#			2) To avoid nearly neutral dynamics, N*m >~ 100.
#		
#			BUT... computation slows dramatically once N > 30,000. The chosen values
#			for N, m, and s are an attempt to optimize these tradeoffs while still 
#			keeping the computation time down to a reasonable duration. They are still
#			not perfect. For example, with s = 0.05 and m = 0.002, N*m = 60; or with
#			s = 0.05, m=0.004, the initial p = 0.16. However, for other parameter 
#			combinations, we satisfy the conditions reasonably well. In a perfect world
#			we would have more time and use larger N and smaller s and m... but so it goes.
#		



rm(list=ls())
#####################
##  Dependencies
source('R/functions-WFSims-Xlinked.R')


######################
##  Run Simulations

# Locally adaptive alleles recessive
	# No deleterious recessive mutations
#	makeFastReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 0, 
#										m.vals = c(0.002, 0.004), m.deltas = NULL,
#										s.vals = c(0.05, 0.1), s.deltas = NULL, 
#										s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
#										r.vals = seq(from = 0, to = 0.5, by = 0.05),
#										newMutant=c("random","random"))


# Additive fitness effects
	# No deleterious recessive mutations
#	makeFastReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 1/2, 
#										m.vals = c(0.0002, 0.0004), m.deltas = NULL,
#										s.vals = c(0.005, 0.01), s.deltas = NULL, 
#										s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
#										r.vals = seq(from = 0, to = 0.5, by = 0.05),
#										newMutant=c("random","random"))

# Locally adaptive alleles dominant
	# No deleterious recessive mutations
#	makeFastReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 1, 
#										m.vals = c(0.0002, 0.0004), m.deltas = NULL,
#										s.vals = c(0.005, 0.01), s.deltas = NULL, 
#										s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
#										r.vals = seq(from = 0, to = 0.5, by = 0.05),
#										newMutant=c("random","random"))


###############################################
# Sims for Fig. 1 etc.

# Locally adaptive alleles recessive
	# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 0, 
									m.vals = c(0.0002), m.deltas = NULL,
									s.vals = c(0.005), s.deltas = NULL, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))

# Additive fitness effects
	# No deleterious recessive mutations
#	makeFigReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 1/2, 
#									m.vals = c(0.0002), m.deltas = NULL,
#									s.vals = c(0.005), s.deltas = NULL, 
#									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
#									r.vals = seq(from = 0, to = 0.5, by = 0.05),
#									newMutant=c("random","random"))

# Locally adaptive alleles dominant
# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = 1000000, N = 500000, h = 1, 
									m.vals = c(0.0002), m.deltas = NULL,
									s.vals = c(0.005), s.deltas = NULL, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))
