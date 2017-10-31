##  Code to run exact simulations for the haploid 2-locus 
##  model of migration-selection balance, and compare
##  the exact genotypic frequencies and LD against the 
##  approximation for the single-locus equilibrium 
##  frequency (as in Kirkpatrick & Barton 2006) 

##  NOTE: To explore the consequences of varying migration rate,
##        r, and recombination rate, r, run this script 
##        interactively, and modify these parameters in the 
##        hapApproxCompare() function. To make a nice pdf of the 
##        resulting figure, run the last couple lines. 
 

rm(list=ls())

## Dependencies
source('./R/functions-figures.R')
source('./R/functions-haploid-exact-eq-sims.R')

##  Run the simulations and display summary plots
res  <-  hapApproxCompare(m=0.05, r=0.5)

##  Make summary plots of [A] and LD
hapApproxComparePlots(data = res)

##  Save plots to .pdf
toPdf(hapApproxComparePlots(), figPath(name='hapApproxCompare.pdf'), width=10, height=5)
embed_fonts(figPath(name='hapApproxCompare.pdf'))