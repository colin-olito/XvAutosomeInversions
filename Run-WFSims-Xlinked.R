################################################################
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



rm(list=ls())
#####################
##  Dependencies
source('R/functions-WFSims-Xlinked.R')


##############################################################################################
##############################################################################################
# Run sims until p*_X corresponding to invasion probability of 0.9997

######################################
##	Small population Size (N = 30k)	##
##	no deleterious mutations       	##
##	m = 0.002, s = 0.05				##
######################################
nReps  <-  1000000
N      <-  30000
m      <-  0.002
s      <-  0.05
# Locally adaptive alleles recessive
	# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 0, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))

# Additive fitness effects
	# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 1/2, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))

# Locally adaptive alleles dominant
# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 1, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))




######################################
##	Large population Size (N = 500k)##
##	no deleterious mutations       	##
##	m = 0.0002, s = 0.005			##
######################################
nReps  <-  1000000
N      <-  500000
m      <-  0.0002
s      <-  0.005
# Locally adaptive alleles recessive
	# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 0, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0,
									newMutant=c("random","random"))

# Additive fitness effects
	# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 1/2, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									newMutant=c("random","random"))

# Locally adaptive alleles dominant
# No deleterious recessive mutations
	makeFigReplicateInvSimsDataXlinked(nReps = nReps, N = N, h = 1, 
									m.vals = m, m.deltas = m,
									s.vals = s, s.deltas = s, 
									s.del.opt = "none", n = 100, u = 1e-5, h.del = 0, 
									r.vals = seq(from = 0, to = 0.5, by = 0.05),
									newMutant=c("random","random"))


##############################################################################################
##############################################################################################
# Run sims until loss of inversion or 4N generations, retain final frequencies

######################################
##	Small population Size (N = 30k)	##
##	Additive fitness LocAdapt fit.	##
##	equal m, s:	m = 0.002, s = 0.05	##
######################################
nReps  <-  50000
N      <-  30000
m      <-  0.002
s      <-  0.05

# No Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "none", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))
# Strong Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "strong", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))
# Lethal Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "lethal", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))


##########################################
##	Large population Size (N = 500k)	##
##	equal m, s:	m = 0.0002, s = 0.005	##
##########################################
nReps  <-  500000
N      <-  500000
m      <-  0.0002
s      <-  0.005

# No Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "none", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))
# Strong Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "strong", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))
# Lethal Deleterious Mutations
makeFigReplicateInvSimsDataXlinked4N(nReps = nReps, N = N, h = 1/2, 
                                     mf = m, mm = m,
                                     sf = s,  sm = s,
                                     r.vals = seq(from = 0, to = 0.5, by = 0.05),
                                     s.del.opt = "lethal", n = 100, u = 1e-5, h.del = 0,
                                     newMutant=c("random","random"))
