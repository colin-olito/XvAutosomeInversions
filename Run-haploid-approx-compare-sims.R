##  Code to run exact simulations for the haploid 2-locus 
##  model of migration-selection balance, and compare
##  the exact genotypic frequencies and LD against the 
##  approximations made in Kirkpatrick & Barton (2006) 

rm(list=ls())

## Dependencies
source('./R/functions-plots.R')
source('./R/functions-haploid-exact-eq-sims.R')

##  Run the simulations and display summary plots
res  <-  hapApproxCompare(m=0.05, r=0.5)

##  Make summary plots of [A] and LD
hapApproxComparePlots(data = res)

##  Save plots to .pdf
toPdf(hapApproxComparePlots(), figPath(name='hapApproxCompare.pdf'), width=10, height=5)
embed_fonts(figPath(name='hapApproxCompare.pdf'))