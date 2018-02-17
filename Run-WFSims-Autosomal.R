##############################################################
#  Wright-Fisher forward simulations of invasion of autosomal
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

rm(list=ls())
#####################
##  Dependencies

source('R/functions-figures.R')
source('R/functions-WFSims-Autosomal.R')


######################
##  Run Simulations

# Locally adaptive alleles recessive
recData  <-  makeReplicateAutoInvSimsData(nReps = 3e+5, N.vals = c(500, 1000), m.vals = c(0.005, 0.01), 
										  s = 0.05, h = 0, r = 0.1, n = 100, u = 1e-5, h.del = 0, 
										  newMutant = 'random') 

# Additive fitness effects
addData  <-  makeReplicateAutoInvSimsData(nReps = 3e+5, N.vals = c(500, 1000), m.vals = c(0.005, 0.01), 
										  s = 0.05, h = 1/2, r = 0.1, n = 100, u = 1e-5, h.del = 0, 
										  newMutant = 'random') 

# Locally adaptive alleles dominant
domData  <-  makeReplicateAutoInvSimsData(nReps = 3e+5, N.vals = c(500, 1000), m.vals = c(0.005, 0.01), 
										  s = 0.05, h = 1, r = 0.1, n = 100, u = 1e-5, h.del = 0, 
										  newMutant = 'random') 
